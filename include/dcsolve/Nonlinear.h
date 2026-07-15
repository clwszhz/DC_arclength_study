#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "dcsolve/Circuit.h"
#include "dcsolve/Lup.h"

namespace dcsolve {

struct NonlinearMnaSystem {
    std::vector<int> node_labels;
    std::vector<std::string> variable_names;
    std::vector<std::string> stamp_log;
};

NonlinearMnaSystem makeNonlinearMnaSystem(const Circuit& circuit);

Vector evaluateResidual(
    const Circuit& circuit, const NonlinearMnaSystem& system, const Vector& x);
Matrix evaluateJacobian(
    const Circuit& circuit, const NonlinearMnaSystem& system, const Vector& x);

void printNonlinearSystem(std::ostream& os, const NonlinearMnaSystem& system);
void printNamedVector(
    std::ostream& os,
    const std::vector<std::string>& names,
    const Vector& values,
    const std::string& title);

} // namespace dcsolve
