#include "dcsolve/Homotopy.h"

#include "dcsolve/Lup.h"
#include "dcsolve/Mna.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dcsolve {
namespace {

double dot(const Vector& a, const Vector& b)
{
    if (a.size() != b.size()) {
        throw std::runtime_error("dot product dimension mismatch");
    }

    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double norm2(const Vector& values)
{
    return std::sqrt(dot(values, values));
}

Vector normalized(Vector values)
{
    const double length = norm2(values);
    if (length == 0.0) {
        throw std::runtime_error("cannot normalize a zero vector");
    }
    for (double& value : values) {
        value /= length;
    }
    return values;
}

Vector negate(Vector values)
{
    for (double& value : values) {
        value = -value;
    }
    return values;
}

Vector makeY(double lambda, const Vector& x)
{
    Vector y;
    y.reserve(x.size() + 1);
    y.push_back(lambda);
    y.insert(y.end(), x.begin(), x.end());
    return y;
}

Vector yToX(const Vector& y)
{
    if (y.empty()) {
        throw std::runtime_error("empty homotopy state");
    }
    return Vector(y.begin() + 1, y.end());
}

Vector addScaled(const Vector& a, double scale, const Vector& b)
{
    if (a.size() != b.size()) {
        throw std::runtime_error("addScaled dimension mismatch");
    }

    Vector result = a;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] += scale * b[i];
    }
    return result;
}

std::size_t findNodeRow(const NonlinearMnaSystem& system, int node)
{
    const auto found =
        std::find(system.node_labels.begin(), system.node_labels.end(), node);
    if (found == system.node_labels.end()) {
        throw std::runtime_error("voltage source node is missing from MNA rows");
    }
    return static_cast<std::size_t>(found - system.node_labels.begin());
}

std::vector<std::size_t> makeHomotopyRowOrder(
    const Circuit& circuit, const NonlinearMnaSystem& system)
{
    const std::size_t n = system.variable_names.size();
    const std::size_t node_count = system.node_labels.size();
    std::vector<std::size_t> rows(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        rows[i] = i;
    }

    for (std::size_t k = 0; k < circuit.voltage_sources.size(); ++k) {
        const auto& source = circuit.voltage_sources[k];
        const std::size_t branch_row = node_count + k;

        int constrained_node = 0;
        if (source.p != 0 && source.n == 0) {
            constrained_node = source.p;
        } else if (source.p == 0 && source.n != 0) {
            constrained_node = source.n;
        } else {
            continue;
        }

        const std::size_t node_row = findNodeRow(system, constrained_node);
        rows[node_row] = branch_row;
        rows[branch_row] = node_row;
    }

    return rows;
}

Vector reorderRows(Vector values, const std::vector<std::size_t>& rows)
{
    if (values.size() != rows.size()) {
        throw std::runtime_error("homotopy row permutation dimension mismatch");
    }

    Vector ordered(values.size(), 0.0);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        ordered[i] = values[rows[i]];
    }
    return ordered;
}

Matrix reorderRows(Matrix values, const std::vector<std::size_t>& rows)
{
    if (values.size() != rows.size()) {
        throw std::runtime_error("homotopy matrix row permutation dimension mismatch");
    }

    Matrix ordered(values.size(), Vector(values.size(), 0.0));
    for (std::size_t i = 0; i < rows.size(); ++i) {
        ordered[i] = values[rows[i]];
    }
    return ordered;
}

Vector evaluateHomotopyResidualBase(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& x)
{
    return reorderRows(
        evaluateResidual(circuit, system, x),
        makeHomotopyRowOrder(circuit, system));
}

Matrix evaluateHomotopyJacobianBase(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& x)
{
    return reorderRows(
        evaluateJacobian(circuit, system, x),
        makeHomotopyRowOrder(circuit, system));
}

bool allFinite(const Vector& values)
{
    for (double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

Matrix evaluateHomotopyJacobianX(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& x,
    double lambda,
    double gleak)
{
    Matrix jacobian = evaluateHomotopyJacobianBase(circuit, system, x);
    const std::size_t n = jacobian.size();

    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            jacobian[row][col] *= lambda;
        }
        jacobian[row][row] += (1.0 - lambda) * gleak;
    }

    return jacobian;
}

Vector evaluateHomotopyJacobianLambda(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    const Vector& x,
    double gleak)
{
    Vector result = evaluateHomotopyResidualBase(circuit, system, x);
    if (result.size() != start.size() || x.size() != start.size()) {
        throw std::runtime_error("homotopy lambda derivative dimension mismatch");
    }

    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] -= gleak * (x[i] - start[i]);
    }
    return result;
}

Vector tangentByFreeColumn(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double lambda,
    const Vector& x,
    double gleak,
    const Vector& previous_tangent);

Vector initialTangent(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double gleak)
{
    if (gleak <= 0.0) {
        throw std::runtime_error("gleak must be positive");
    }

    Vector previous(start.size() + 1, 0.0);
    previous[0] = 1.0;
    try {
        return tangentByFreeColumn(
            circuit, system, start, 0.0, start, gleak, previous);
    } catch (const std::exception&) {
        const Vector d_h_d_lambda =
            evaluateHomotopyJacobianLambda(circuit, system, start, start, gleak);
        Vector tangent(d_h_d_lambda.size() + 1, 0.0);
        tangent[0] = 1.0;
        for (std::size_t i = 0; i < d_h_d_lambda.size(); ++i) {
            tangent[i + 1] = -d_h_d_lambda[i] / gleak;
        }
        return normalized(std::move(tangent));
    }
}

Matrix makeAugmentedJacobian(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double lambda,
    const Vector& x,
    double gleak,
    const Vector& tangent)
{
    const Matrix d_h_d_x =
        evaluateHomotopyJacobianX(circuit, system, x, lambda, gleak);
    const Vector d_h_d_lambda =
        evaluateHomotopyJacobianLambda(circuit, system, start, x, gleak);

    const std::size_t n = x.size();
    Matrix augmented(n + 1, Vector(n + 1, 0.0));
    for (std::size_t row = 0; row < n; ++row) {
        augmented[row][0] = d_h_d_lambda[row];
        for (std::size_t col = 0; col < n; ++col) {
            augmented[row][col + 1] = d_h_d_x[row][col];
        }
    }
    augmented[n] = tangent;
    return augmented;
}

Vector tangentByFreeColumn(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double lambda,
    const Vector& x,
    double gleak,
    const Vector& previous_tangent)
{
    const Matrix d_h_d_x =
        evaluateHomotopyJacobianX(circuit, system, x, lambda, gleak);
    const Vector d_h_d_lambda =
        evaluateHomotopyJacobianLambda(circuit, system, start, x, gleak);

    const std::size_t n = x.size();
    const std::size_t free_col = n;
    Matrix coefficient(n, Vector(n, 0.0));
    Vector rhs(n, 0.0);

    for (std::size_t row = 0; row < n; ++row) {
        std::size_t out_col = 0;
        for (std::size_t col = 0; col < n + 1; ++col) {
            const double value =
                col == 0 ? d_h_d_lambda[row] : d_h_d_x[row][col - 1];
            if (col == free_col) {
                rhs[row] = -value;
            } else {
                coefficient[row][out_col] = value;
                ++out_col;
            }
        }
    }

    const Vector solved = solveLup(std::move(coefficient), rhs);
    Vector tangent(n + 1, 0.0);
    std::size_t in_col = 0;
    for (std::size_t col = 0; col < n + 1; ++col) {
        if (col == free_col) {
            tangent[col] = 1.0;
        } else {
            tangent[col] = solved[in_col];
            ++in_col;
        }
    }

    tangent = normalized(std::move(tangent));
    if (dot(tangent, previous_tangent) < 0.0) {
        tangent = negate(std::move(tangent));
    }
    return tangent;
}

struct CorrectorResult {
    Vector y;
    std::size_t iterations = 0;
    double norm = 0.0;
    bool converged = false;
};

double pseudoResidualNorm(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    const Vector& y_predicted,
    const Vector& y,
    const Vector& tangent,
    const HomotopyOptions& options)
{
    Vector residual = evaluateHomotopy(
        circuit, system, start, y[0], yToX(y), options.gleak);
    residual.push_back(dot(tangent, addScaled(y, -1.0, y_predicted)));
    return maxAbs(residual);
}

CorrectorResult correctPseudoArclength(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    const Vector& y_predicted,
    Vector y,
    const Vector& tangent,
    const HomotopyOptions& options)
{
    const std::size_t n = start.size();
    if (y.size() != n + 1 || y_predicted.size() != n + 1
        || tangent.size() != n + 1) {
        throw std::runtime_error("pseudo-arclength corrector dimension mismatch");
    }

    for (std::size_t iter = 0; iter <= options.max_newton_iterations; ++iter) {
        const double lambda = y[0];
        const Vector x = yToX(y);
        Vector residual = evaluateHomotopy(
            circuit, system, start, lambda, x, options.gleak);
        residual.push_back(dot(tangent, addScaled(y, -1.0, y_predicted)));

        const double residual_norm = maxAbs(residual);
        if (residual_norm <= options.tolerance) {
            return {y, iter, residual_norm, true};
        }

        if (iter == options.max_newton_iterations) {
            break;
        }

        Matrix augmented = makeAugmentedJacobian(
            circuit, system, start, lambda, x, options.gleak, tangent);
        for (double& value : residual) {
            value = -value;
        }

        const Vector delta = solveLup(std::move(augmented), residual);
        y = addScaled(y, 1.0, delta);
    }

    return {y,
            options.max_newton_iterations,
            pseudoResidualNorm(
                circuit, system, start, y_predicted, y, tangent, options),
            false};
}

Vector tangentAt(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    const Vector& y,
    const Vector& previous_tangent,
    double gleak)
{
    const double lambda = y[0];
    const Vector x = yToX(y);
    try {
        return tangentByFreeColumn(
            circuit, system, start, lambda, x, gleak, previous_tangent);
    } catch (const std::exception&) {
        // Fall back to an augmented solve if the chosen free column is singular.
    }

    Matrix augmented = makeAugmentedJacobian(
        circuit, system, start, lambda, x, gleak, previous_tangent);

    Vector rhs(x.size() + 1, 0.0);
    rhs.back() = 1.0;

    Vector tangent = solveLup(std::move(augmented), rhs);
    tangent = normalized(std::move(tangent));
    if (dot(tangent, previous_tangent) < 0.0) {
        tangent = negate(std::move(tangent));
    }
    return tangent;
}

struct FinalNewtonResult {
    Vector x;
    double norm = 0.0;
    std::size_t iterations = 0;
    bool converged = false;
};

FinalNewtonResult solveAtLambdaOne(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    Vector x,
    const HomotopyOptions& options)
{
    for (std::size_t iter = 0; iter <= options.max_final_newton_iterations; ++iter) {
        Vector residual = evaluateResidual(circuit, system, x);
        const double residual_norm = maxAbs(residual);
        if (residual_norm <= options.tolerance) {
            return {x, residual_norm, iter, true};
        }

        if (iter == options.max_final_newton_iterations) {
            return {x, residual_norm, iter, false};
        }

        Matrix jacobian = evaluateJacobian(circuit, system, x);
        for (double& value : residual) {
            value = -value;
        }
        Vector delta;
        try {
            delta = solveLup(std::move(jacobian), residual);
        } catch (const std::exception&) {
            return {x, residual_norm, iter, false};
        }
        x = addScaled(x, 1.0, delta);
    }

    return {x, 0.0, 0, false};
}

bool sameVariables(
    const NonlinearMnaSystem& nonlinear_system, const MnaSystem& linear_system)
{
    return nonlinear_system.variable_names == linear_system.variable_names;
}

} // namespace

Vector makeDefaultHomotopyStart(
    const Circuit& circuit, const NonlinearMnaSystem& system)
{
    if (!circuit.initial_guess.empty()) {
        if (circuit.initial_guess.size() != system.variable_names.size()) {
            throw std::runtime_error(
                ".START value count must match the nonlinear unknown count");
        }
        return circuit.initial_guess;
    }

    Circuit linear_part = circuit;
    linear_part.bjts.clear();
    linear_part.bjt_models.clear();

    try {
        const MnaSystem linear_system = assembleMna(linear_part);
        if (sameVariables(system, linear_system)) {
            return solveLup(linear_system.a, linear_system.b);
        }
    } catch (const std::exception&) {
        // A zero start still gives H(0, x)=0. The linear start is only a helper.
    }

    return Vector(system.variable_names.size(), 0.0);
}

Vector evaluateHomotopy(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double lambda,
    const Vector& x,
    double gleak)
{
    if (start.size() != x.size()) {
        throw std::runtime_error("homotopy dimension mismatch");
    }

    Vector h = evaluateHomotopyResidualBase(circuit, system, x);
    for (std::size_t i = 0; i < h.size(); ++i) {
        h[i] = (1.0 - lambda) * gleak * (x[i] - start[i]) + lambda * h[i];
    }
    return h;
}

HomotopyResult solvePseudoArclength(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    Vector start,
    const HomotopyOptions& options)
{
    if (options.step_size <= 0.0) {
        throw std::runtime_error("step_size must be positive");
    }
    if (options.gleak <= 0.0) {
        throw std::runtime_error("gleak must be positive");
    }
    if (start.empty()) {
        start = makeDefaultHomotopyStart(circuit, system);
    }
    if (start.size() != system.variable_names.size()) {
        throw std::runtime_error("homotopy start dimension mismatch");
    }

    HomotopyResult result;
    result.steps.push_back({0,
                            0.0,
                            start,
                            maxAbs(evaluateHomotopy(circuit,
                                                    system,
                                                    start,
                                                    0.0,
                                                    start,
                                                    options.gleak)),
                            maxAbs(evaluateResidual(circuit, system, start)),
                            0});

    Vector y = makeY(0.0, start);
    Vector tangent = initialTangent(circuit, system, start, options.gleak);

    double arc_step = options.step_size;
    constexpr double min_arc_step = 1e-8;
    for (std::size_t step = 1; step <= options.max_steps; ++step) {
        CorrectorResult corrected;
        Vector y_predicted;
        bool accepted = false;

        while (arc_step >= min_arc_step) {
            y_predicted = addScaled(y, arc_step, tangent);
            try {
                corrected = correctPseudoArclength(
                    circuit, system, start, y_predicted, y_predicted, tangent, options);
            } catch (const std::exception&) {
                corrected.converged = false;
            }

            if (corrected.converged && allFinite(corrected.y)) {
                accepted = true;
                break;
            }

            arc_step *= 0.5;
        }

        if (!accepted) {
            break;
        }

        y = corrected.y;
        const Vector x = yToX(y);
        result.steps.push_back({step,
                                y[0],
                                x,
                                corrected.norm,
                                maxAbs(evaluateResidual(circuit, system, x)),
                                corrected.iterations});

        if (y[0] >= 1.0 - options.tolerance) {
            const FinalNewtonResult final =
                solveAtLambdaOne(circuit, system, x, options);
            result.converged = final.converged;
            result.lambda = 1.0;
            result.x = final.x;
            result.steps.push_back({step + 1,
                                    1.0,
                                    final.x,
                                    maxAbs(evaluateHomotopy(circuit,
                                                            system,
                                                            start,
                                                            1.0,
                                                            final.x,
                                                            options.gleak)),
                                    final.norm,
                                    final.iterations});
            return result;
        }

        try {
            tangent = tangentAt(
                circuit, system, start, y, tangent, options.gleak);
        } catch (const std::exception&) {
            arc_step *= 0.5;
        }

        if (corrected.iterations <= 3) {
            arc_step = std::min(options.step_size, arc_step * 1.25);
        }
    }

    result.converged = false;
    result.lambda = y[0];
    result.x = yToX(y);
    return result;
}

} // namespace dcsolve
