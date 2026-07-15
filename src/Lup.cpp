#include "dcsolve/Lup.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dcsolve {
namespace {

void checkSquare(const Matrix& a)
{
    const std::size_t n = a.size();
    for (const auto& row : a) {
        if (row.size() != n) {
            throw std::runtime_error("LUP requires a square matrix");
        }
    }
}

} // namespace

LupDecomposition factorizeLup(Matrix a, double tolerance)
{
    checkSquare(a);

    const std::size_t n = a.size();
    LupDecomposition result;
    result.lu = std::move(a);
    result.pivots.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        result.pivots[i] = i;
    }

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot_row = k;
        double pivot_abs = std::abs(result.lu[k][k]);

        for (std::size_t row = k + 1; row < n; ++row) {
            const double candidate = std::abs(result.lu[row][k]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }

        if (pivot_abs <= tolerance) {
            throw std::runtime_error("matrix is singular or nearly singular at column "
                                     + std::to_string(k));
        }

        if (pivot_row != k) {
            std::swap(result.lu[pivot_row], result.lu[k]);
            std::swap(result.pivots[pivot_row], result.pivots[k]);
            ++result.swap_count;
        }

        for (std::size_t row = k + 1; row < n; ++row) {
            result.lu[row][k] /= result.lu[k][k];
            for (std::size_t col = k + 1; col < n; ++col) {
                result.lu[row][col] -= result.lu[row][k] * result.lu[k][col];
            }
        }
    }

    return result;
}

Vector solveLup(const LupDecomposition& decomp, const Vector& b)
{
    const std::size_t n = decomp.lu.size();
    if (b.size() != n || decomp.pivots.size() != n) {
        throw std::runtime_error("LUP solve dimension mismatch");
    }

    Vector y(n, 0.0);
    Vector x(n, 0.0);

    // Forward substitution: L*y = P*b. L has an implicit unit diagonal.
    for (std::size_t i = 0; i < n; ++i) {
        double sum = b[decomp.pivots[i]];
        for (std::size_t j = 0; j < i; ++j) {
            sum -= decomp.lu[i][j] * y[j];
        }
        y[i] = sum;
    }

    // Back substitution: U*x = y.
    for (std::size_t reverse_i = 0; reverse_i < n; ++reverse_i) {
        const std::size_t i = n - 1 - reverse_i;
        double sum = y[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= decomp.lu[i][j] * x[j];
        }
        x[i] = sum / decomp.lu[i][i];
    }

    return x;
}

Vector solveLup(Matrix a, const Vector& b, double tolerance)
{
    return solveLup(factorizeLup(std::move(a), tolerance), b);
}

} // namespace dcsolve
