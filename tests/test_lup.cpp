#include "dcsolve/Circuit.h"
#include "dcsolve/Homotopy.h"
#include "dcsolve/Lup.h"
#include "dcsolve/Mna.h"
#include "dcsolve/Nonlinear.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void requireNear(double actual, double expected, double tolerance, const char* name)
{
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(std::string(name) + " expected "
                                 + std::to_string(expected) + " but got "
                                 + std::to_string(actual));
    }
}

void testClassic3x3()
{
    dcsolve::Matrix a = {
        {2.0, 1.0, -1.0},
        {-3.0, -1.0, 2.0},
        {-2.0, 1.0, 2.0},
    };
    dcsolve::Vector b = {8.0, -11.0, -3.0};

    const auto x = dcsolve::solveLup(a, b);
    requireNear(x[0], 2.0, 1e-10, "x0");
    requireNear(x[1], 3.0, 1e-10, "x1");
    requireNear(x[2], -1.0, 1e-10, "x2");
}

void testPivoting()
{
    dcsolve::Matrix a = {
        {0.0, 2.0},
        {1.0, 1.0},
    };
    dcsolve::Vector b = {4.0, 3.0};

    const auto x = dcsolve::solveLup(a, b);
    requireNear(x[0], 1.0, 1e-10, "pivot x0");
    requireNear(x[1], 2.0, 1e-10, "pivot x1");
}

void testDividerMna()
{
    dcsolve::Circuit circuit;
    circuit.has_op = true;
    circuit.voltage_sources.push_back({"V1", 1, 0, 10.0});
    circuit.resistors.push_back({"R1", 1, 2, 1000.0});
    circuit.resistors.push_back({"R2", 2, 0, 1000.0});

    const auto system = dcsolve::assembleMna(circuit);
    const auto x = dcsolve::solveLup(system.a, system.b);

    requireNear(x[0], 10.0, 1e-10, "divider V(1)");
    requireNear(x[1], 5.0, 1e-10, "divider V(2)");
    requireNear(x[2], -0.005, 1e-10, "divider I(V1)");

    const auto r = dcsolve::residual(system.a, x, system.b);
    requireNear(dcsolve::maxAbs(r), 0.0, 1e-12, "divider residual");
}

dcsolve::Circuit makeDividerCircuit()
{
    dcsolve::Circuit circuit;
    circuit.has_op = true;
    circuit.voltage_sources.push_back({"V1", 1, 0, 10.0});
    circuit.resistors.push_back({"R1", 1, 2, 1000.0});
    circuit.resistors.push_back({"R2", 2, 0, 1000.0});
    return circuit;
}

void testNonlinearDividerHomotopy()
{
    const auto circuit = makeDividerCircuit();
    const auto system = dcsolve::makeNonlinearMnaSystem(circuit);

    dcsolve::HomotopyOptions options;
    options.step_size = 0.5;
    options.max_steps = 4;

    const auto result =
        dcsolve::solvePseudoArclength(circuit, system, {}, options);
    if (!result.converged) {
        throw std::runtime_error("nonlinear divider homotopy did not converge");
    }

    requireNear(result.x[0], 10.0, 1e-9, "nonlinear divider V(1)");
    requireNear(result.x[1], 5.0, 1e-9, "nonlinear divider V(2)");
    requireNear(result.x[2], -0.005, 1e-9, "nonlinear divider I(V1)");

    const auto f = dcsolve::evaluateResidual(circuit, system, result.x);
    requireNear(dcsolve::maxAbs(f), 0.0, 1e-10, "nonlinear divider F");
}

void testBjtJacobianAgainstFiniteDifference()
{
    dcsolve::Circuit circuit;
    circuit.has_op = true;
    circuit.bjt_models.push_back({"QN", dcsolve::BjtType::Npn, 1e-16, 0.99, 0.5, 38.78});
    circuit.bjts.push_back({"Q1", 1, 2, 3, "QN"});

    const auto system = dcsolve::makeNonlinearMnaSystem(circuit);
    const dcsolve::Vector x = {1.0, 0.7, 0.0};
    const auto analytic = dcsolve::evaluateJacobian(circuit, system, x);

    const double eps = 1e-6;
    for (std::size_t col = 0; col < x.size(); ++col) {
        dcsolve::Vector xp = x;
        dcsolve::Vector xm = x;
        xp[col] += eps;
        xm[col] -= eps;

        const auto fp = dcsolve::evaluateResidual(circuit, system, xp);
        const auto fm = dcsolve::evaluateResidual(circuit, system, xm);
        for (std::size_t row = 0; row < x.size(); ++row) {
            const double numeric = (fp[row] - fm[row]) / (2.0 * eps);
            requireNear(
                analytic[row][col], numeric, 1e-7, "BJT Jacobian finite diff");
        }
    }
}

} // namespace

int main()
{
    try {
        testClassic3x3();
        testPivoting();
        testDividerMna();
        testNonlinearDividerHomotopy();
        testBjtJacobianAgainstFiniteDifference();
    } catch (const std::exception& e) {
        std::cerr << "test failed: " << e.what() << '\n';
        return 1;
    }

    std::cout << "all tests passed\n";
    return 0;
}
