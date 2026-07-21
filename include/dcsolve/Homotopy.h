#pragma once

#include <cstddef>
#include <vector>

#include "dcsolve/Circuit.h"
#include "dcsolve/Nonlinear.h"

namespace dcsolve {

struct HomotopyOptions {
    double step_size = 0.1;
    double gleak = 1e-3;
    double tolerance = 1e-9;
    std::size_t max_steps = 3000;
    std::size_t max_newton_iterations = 12;
    std::size_t max_final_newton_iterations = 20;
};

struct HomotopyStep {
    std::size_t step = 0;
    double lambda = 0.0;
    Vector x;
    double homotopy_norm = 0.0;
    double circuit_norm = 0.0;
    std::size_t newton_iterations = 0;
};

struct HomotopySolution {
    std::size_t step = 0;
    Vector x;
    double residual_norm = 0.0;
    std::size_t newton_iterations = 0;
};

struct HomotopyResult {
    bool converged = false;
    double lambda = 0.0;
    Vector x;
    std::vector<HomotopyStep> steps;
    std::vector<HomotopySolution> solutions;
};

Vector makeDefaultHomotopyStart(
    const Circuit& circuit, const NonlinearMnaSystem& system);

Vector evaluateHomotopy(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    const Vector& start,
    double lambda,
    const Vector& x,
    double gleak);

HomotopyResult solvePseudoArclength(
    const Circuit& circuit,
    const NonlinearMnaSystem& system,
    Vector start,
    const HomotopyOptions& options = {});

} // namespace dcsolve
