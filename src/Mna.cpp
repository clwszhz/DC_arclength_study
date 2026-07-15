#include "dcsolve/Mna.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace dcsolve {
namespace {

std::string nodeText(int node)
{
    return node == 0 ? "0(GND)" : std::to_string(node);
}

std::string formatDouble(double value)
{
    std::ostringstream out;
    out << std::scientific << std::setprecision(6) << value;
    return out.str();
}

std::unordered_map<int, std::size_t> makeNodeIndex(const std::vector<int>& labels)
{
    std::unordered_map<int, std::size_t> index;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        index[labels[i]] = i;
    }
    return index;
}

void addIfNode(
    Matrix& a,
    const std::unordered_map<int, std::size_t>& node_index,
    int row_node,
    std::size_t col,
    double value)
{
    if (row_node == 0) {
        return;
    }
    a.at(node_index.at(row_node)).at(col) += value;
}

void addRhsIfNode(
    Vector& b,
    const std::unordered_map<int, std::size_t>& node_index,
    int node,
    double value)
{
    if (node == 0) {
        return;
    }
    b.at(node_index.at(node)) += value;
}

} // namespace

MnaSystem assembleMna(const Circuit& circuit)
{
    if (!circuit.bjts.empty()) {
        throw std::runtime_error(
            "linear MNA does not support BJT elements; use homotopy_solve");
    }

    MnaSystem system;
    system.node_labels = circuit.nodeLabels();

    const auto node_index = makeNodeIndex(system.node_labels);
    const std::size_t node_count = system.node_labels.size();
    const std::size_t voltage_count = circuit.voltage_sources.size();
    const std::size_t n = node_count + voltage_count;

    system.variable_names.reserve(n);
    for (int label : system.node_labels) {
        system.variable_names.push_back("V(" + std::to_string(label) + ")");
    }
    for (const auto& source : circuit.voltage_sources) {
        system.variable_names.push_back("I(" + source.name + ")");
    }

    system.a.assign(n, Vector(n, 0.0));
    system.b.assign(n, 0.0);

    for (const auto& r : circuit.resistors) {
        const double g = 1.0 / r.resistance;

        if (r.a != 0) {
            const std::size_t ia = node_index.at(r.a);
            system.a[ia][ia] += g;
        }
        if (r.b != 0) {
            const std::size_t ib = node_index.at(r.b);
            system.a[ib][ib] += g;
        }
        if (r.a != 0 && r.b != 0) {
            const std::size_t ia = node_index.at(r.a);
            const std::size_t ib = node_index.at(r.b);
            system.a[ia][ib] -= g;
            system.a[ib][ia] -= g;
        }

        system.stamp_log.push_back(
            r.name + ": resistor between node " + nodeText(r.a) + " and "
            + nodeText(r.b) + ", g=1/R=" + formatDouble(g)
            + " stamps +g on diagonals and -g on the two off-diagonals");
    }

    for (const auto& source : circuit.current_sources) {
        // Positive source current flows from p to n. In KCL form A*x=b,
        // it leaves node p and enters node n.
        addRhsIfNode(system.b, node_index, source.p, -source.current);
        addRhsIfNode(system.b, node_index, source.n, source.current);

        system.stamp_log.push_back(
            source.name + ": current source " + formatDouble(source.current)
            + " A from node " + nodeText(source.p) + " to "
            + nodeText(source.n)
            + " stamps b[p]-=I and b[n]+=I");
    }

    for (std::size_t k = 0; k < circuit.voltage_sources.size(); ++k) {
        const auto& source = circuit.voltage_sources[k];
        const std::size_t col = node_count + k;
        const std::size_t row = node_count + k;

        addIfNode(system.a, node_index, source.p, col, 1.0);
        addIfNode(system.a, node_index, source.n, col, -1.0);

        if (source.p != 0) {
            system.a[row][node_index.at(source.p)] += 1.0;
        }
        if (source.n != 0) {
            system.a[row][node_index.at(source.n)] -= 1.0;
        }
        system.b[row] += source.voltage;

        system.stamp_log.push_back(
            source.name + ": voltage source enforces V("
            + std::to_string(source.p) + ")-V(" + std::to_string(source.n)
            + ")=" + formatDouble(source.voltage)
            + " and introduces unknown I(" + source.name + ")");
    }

    if (system.a.empty()) {
        throw std::runtime_error("empty MNA system: add at least one non-ground node");
    }

    return system;
}

Vector residual(const Matrix& a, const Vector& x, const Vector& b)
{
    if (a.size() != b.size() || x.size() != b.size()) {
        throw std::runtime_error("residual dimension mismatch");
    }

    Vector r(b.size(), 0.0);
    for (std::size_t row = 0; row < a.size(); ++row) {
        if (a[row].size() != x.size()) {
            throw std::runtime_error("residual matrix is not rectangular");
        }
        for (std::size_t col = 0; col < x.size(); ++col) {
            r[row] += a[row][col] * x[col];
        }
        r[row] -= b[row];
    }
    return r;
}

double maxAbs(const Vector& values)
{
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

void printMnaSystem(std::ostream& os, const MnaSystem& system)
{
    os << "Variable order x:\n";
    for (std::size_t i = 0; i < system.variable_names.size(); ++i) {
        os << "  x[" << i << "] = " << system.variable_names[i] << '\n';
    }

    os << "\nStamp log:\n";
    for (const auto& line : system.stamp_log) {
        os << "  - " << line << '\n';
    }

    os << "\nMNA system A*x=b:\n";
    os << std::scientific << std::setprecision(6);
    for (std::size_t row = 0; row < system.a.size(); ++row) {
        os << "  [";
        for (std::size_t col = 0; col < system.a[row].size(); ++col) {
            os << std::setw(14) << system.a[row][col];
            if (col + 1 != system.a[row].size()) {
                os << ' ';
            }
        }
        os << " ]  [x] = " << std::setw(14) << system.b[row] << '\n';
    }
}

void printSolution(std::ostream& os, const MnaSystem& system, const Vector& x)
{
    os << "\nSolution:\n";
    os << std::scientific << std::setprecision(9);
    for (std::size_t i = 0; i < x.size(); ++i) {
        const bool is_voltage = system.variable_names[i].rfind("V(", 0) == 0;
        os << "  " << std::setw(10) << system.variable_names[i] << " = "
           << std::setw(16) << x[i] << (is_voltage ? " V" : " A") << '\n';
    }
}

void printResidual(std::ostream& os, const Vector& r)
{
    os << "\nResidual r=A*x-b:\n";
    os << std::scientific << std::setprecision(9);
    for (std::size_t i = 0; i < r.size(); ++i) {
        os << "  r[" << i << "] = " << r[i] << '\n';
    }
    os << "  max|r| = " << maxAbs(r) << '\n';
}

} // namespace dcsolve
