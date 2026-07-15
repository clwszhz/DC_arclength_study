#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "dcsolve/Circuit.h"
#include "dcsolve/Lup.h"

namespace dcsolve {

struct MnaSystem {
    std::vector<int> node_labels;
    std::vector<std::string> variable_names;
    Matrix a;
    Vector b;
    std::vector<std::string> stamp_log;
};

MnaSystem assembleMna(const Circuit& circuit);

Vector residual(const Matrix& a, const Vector& x, const Vector& b);
double maxAbs(const Vector& values);

void printMnaSystem(std::ostream& os, const MnaSystem& system);
void printSolution(std::ostream& os, const MnaSystem& system, const Vector& x);
void printResidual(std::ostream& os, const Vector& r);

} // namespace dcsolve
