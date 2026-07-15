#include "dcsolve/Circuit.h"
#include "dcsolve/Lup.h"
#include "dcsolve/Mna.h"

#include <exception>
#include <iostream>
#include <string>

namespace {

void printUsage(const char* program)
{
    std::cerr << "usage: " << program << " <netlist.cir>\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        printUsage(argv[0]);
        return 2;
    }

    try {
        const std::string path = argv[1];
        const auto circuit = dcsolve::parseNetlistFile(path);
        const auto system = dcsolve::assembleMna(circuit);

        dcsolve::printMnaSystem(std::cout, system);

        const auto solution = dcsolve::solveLup(system.a, system.b);
        dcsolve::printSolution(std::cout, system, solution);

        const auto r = dcsolve::residual(system.a, solution, system.b);
        dcsolve::printResidual(std::cout, r);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
