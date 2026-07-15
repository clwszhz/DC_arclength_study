#pragma once

#include <cstddef>
#include <vector>

namespace dcsolve {

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

struct LupDecomposition {
    Matrix lu;
    std::vector<std::size_t> pivots;
    int swap_count = 0;
};

LupDecomposition factorizeLup(Matrix a, double tolerance = 1e-12);
Vector solveLup(const LupDecomposition& decomp, const Vector& b);
Vector solveLup(Matrix a, const Vector& b, double tolerance = 1e-12);

} // namespace dcsolve
