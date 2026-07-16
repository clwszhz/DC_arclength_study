#include "dcsolve/Circuit.h"
#include "dcsolve/Homotopy.h"
#include "dcsolve/Lup.h"
#include "dcsolve/Mna.h"
#include "dcsolve/Nonlinear.h"

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

void printUsage(const char* program)
{
    std::cerr << "usage: " << program
              << " <netlist.cir> [step_size] [max_steps] [gleak]\n";
    std::cerr << "  R/I/V-only netlists use direct linear MNA.\n";
    std::cerr << "  Netlists with BJT elements use nonlinear homotopy.\n";
}

double parseDoubleArg(const char* text, const char* name)
{
    try {
        return std::stod(text);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string(name) + " must be a number");
    }
}

std::size_t parseSizeArg(const char* text, const char* name)
{
    try {
        return static_cast<std::size_t>(std::stoul(text));
    } catch (const std::exception&) {
        throw std::runtime_error(std::string(name) + " must be an integer");
    }
}

bool hasNonlinearDevices(const dcsolve::Circuit& circuit)
{
    return !circuit.bjts.empty();
}

void printProgress(const dcsolve::HomotopyResult& result)
{
    std::cout << "\nHomotopy progress:\n";
    std::cout << "  " << std::setw(5) << "step" << ' ' << std::setw(12)
              << "lambda" << ' ' << std::setw(14) << "max|H|" << ' '
              << std::setw(14) << "max|F|" << ' ' << std::setw(8)
              << "newton" << '\n';

    std::cout << std::scientific << std::setprecision(6);
    for (const auto& step : result.steps) {
        std::cout << "  " << std::setw(5) << step.step << ' ' << std::setw(12)
                  << step.lambda << ' ' << std::setw(14)
                  << step.homotopy_norm << ' ' << std::setw(14)
                  << step.circuit_norm << ' ' << std::setw(8)
                  << step.newton_iterations << '\n';
    }
}

int solveLinear(const dcsolve::Circuit& circuit)
{
    std::cout << "Analysis mode: direct linear MNA\n";
    const auto system = dcsolve::assembleMna(circuit);

    dcsolve::printMnaSystem(std::cout, system);

    const auto solution = dcsolve::solveLup(system.a, system.b);
    dcsolve::printSolution(std::cout, system, solution);

    const auto r = dcsolve::residual(system.a, solution, system.b);
    dcsolve::printResidual(std::cout, r);
    return 0;
}

int solveNonlinear(const dcsolve::Circuit& circuit, const dcsolve::HomotopyOptions& options)
{
    std::cout << "Analysis mode: nonlinear homotopy\n";
    const auto system = dcsolve::makeNonlinearMnaSystem(circuit);
    dcsolve::printNonlinearSystem(std::cout, system);

    const auto start = dcsolve::makeDefaultHomotopyStart(circuit, system);
    dcsolve::printNamedVector(
        std::cout, system.variable_names, start, "Homotopy start a");

    const auto result =
        dcsolve::solvePseudoArclength(circuit, system, start, options);
    printProgress(result);

    dcsolve::printNamedVector(
        std::cout, system.variable_names, result.x, "Final solution");
    const auto final_residual =
        dcsolve::evaluateResidual(circuit, system, result.x);
    dcsolve::printNamedVector(
        std::cout, system.variable_names, final_residual, "Final F(x)");
    std::cout << "  max|F| = " << dcsolve::maxAbs(final_residual) << '\n';

    if (!result.converged) {
        std::cerr << "warning: homotopy did not converge to lambda=1\n";
        return 1;
    }

    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 5) {
        printUsage(argv[0]);
        return 2;
    }

    try {
        dcsolve::HomotopyOptions options;
        if (argc >= 3) {
            options.step_size = parseDoubleArg(argv[2], "step_size");
        }
        if (argc >= 4) {
            options.max_steps = parseSizeArg(argv[3], "max_steps");
        }
        if (argc >= 5) {
            options.gleak = parseDoubleArg(argv[4], "gleak");
        }

        const std::string path = argv[1];
        const auto circuit = dcsolve::parseNetlistFile(path);

        if (hasNonlinearDevices(circuit)) {
            return solveNonlinear(circuit, options);
        }
        return solveLinear(circuit);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
