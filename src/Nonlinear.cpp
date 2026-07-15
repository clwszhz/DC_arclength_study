#include "dcsolve/Nonlinear.h"

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

double nodeVoltage(
    const Vector& x,
    const std::unordered_map<int, std::size_t>& node_index,
    int node)
{
    if (node == 0) {
        return 0.0;
    }
    return x.at(node_index.at(node));
}

void addResidualIfNode(
    Vector& f,
    const std::unordered_map<int, std::size_t>& node_index,
    int node,
    double value)
{
    if (node == 0) {
        return;
    }
    f.at(node_index.at(node)) += value;
}

void addJacobianNodeColumnIfNode(
    Matrix& jacobian,
    const std::unordered_map<int, std::size_t>& node_index,
    int row_node,
    int col_node,
    double value)
{
    if (row_node == 0 || col_node == 0) {
        return;
    }
    jacobian.at(node_index.at(row_node)).at(node_index.at(col_node)) += value;
}

void addJacobianExplicitColumnIfNode(
    Matrix& jacobian,
    const std::unordered_map<int, std::size_t>& node_index,
    int row_node,
    std::size_t col,
    double value)
{
    if (row_node == 0) {
        return;
    }
    jacobian.at(node_index.at(row_node)).at(col) += value;
}

const BjtModel& findBjtModel(const Circuit& circuit, const std::string& model_name)
{
    for (const auto& model : circuit.bjt_models) {
        if (model.name == model_name) {
            return model;
        }
    }
    throw std::runtime_error("missing BJT model: " + model_name);
}

double limitedExp(double value)
{
    return std::exp(std::clamp(value, -80.0, 80.0));
}

struct BjtStamp {
    double ic = 0.0;
    double ib = 0.0;
    double ie = 0.0;

    double dic_dvc = 0.0;
    double dic_dvb = 0.0;
    double dic_dve = 0.0;
    double dib_dvc = 0.0;
    double dib_dvb = 0.0;
    double dib_dve = 0.0;
    double die_dvc = 0.0;
    double die_dvb = 0.0;
    double die_dve = 0.0;
};

BjtStamp evaluateBjtStamp(
    const BjtModel& model, double vc, double vb, double ve)
{
    const double is = model.saturation_current;
    const double af = model.alpha_forward;
    const double ar = model.alpha_reverse;
    const double n = model.exponential_coefficient;

    const double ebe = limitedExp(n * (vb - ve));
    const double ebc = limitedExp(n * (vb - vc));

    const double fe = is / af * (ebe - 1.0);
    const double fc = is / ar * (ebc - 1.0);

    const double dfe_dvb = is / af * n * ebe;
    const double dfe_dve = -dfe_dvb;
    const double dfc_dvb = is / ar * n * ebc;
    const double dfc_dvc = -dfc_dvb;

    BjtStamp stamp;
    stamp.ie = fe - ar * fc;
    stamp.ic = fc - af * fe;
    stamp.ib = -stamp.ic - stamp.ie;

    stamp.die_dvc = -ar * dfc_dvc;
    stamp.die_dvb = dfe_dvb - ar * dfc_dvb;
    stamp.die_dve = dfe_dve;

    stamp.dic_dvc = dfc_dvc;
    stamp.dic_dvb = dfc_dvb - af * dfe_dvb;
    stamp.dic_dve = -af * dfe_dve;

    stamp.dib_dvc = -stamp.dic_dvc - stamp.die_dvc;
    stamp.dib_dvb = -stamp.dic_dvb - stamp.die_dvb;
    stamp.dib_dve = -stamp.dic_dve - stamp.die_dve;
    return stamp;
}

} // namespace

NonlinearMnaSystem makeNonlinearMnaSystem(const Circuit& circuit)
{
    NonlinearMnaSystem system;
    system.node_labels = circuit.nodeLabels();

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

    for (const auto& r : circuit.resistors) {
        const double g = 1.0 / r.resistance;
        system.stamp_log.push_back(
            r.name + ": nonlinear F/J resistor stamp between node "
            + nodeText(r.a) + " and " + nodeText(r.b)
            + ", current is g*(Va-Vb), g=" + formatDouble(g));
    }
    for (const auto& source : circuit.current_sources) {
        system.stamp_log.push_back(
            source.name + ": current source from node " + nodeText(source.p)
            + " to " + nodeText(source.n)
            + " adds +I to F[p] and -I to F[n]");
    }
    for (const auto& source : circuit.voltage_sources) {
        system.stamp_log.push_back(
            source.name + ": voltage source adds branch current and constraint V(p)-V(n)-V=0");
    }
    for (const auto& q : circuit.bjts) {
        system.stamp_log.push_back(
            q.name + ": NPN Ebers-Moll BJT, nodes C=" + nodeText(q.collector)
            + ", B=" + nodeText(q.base) + ", E=" + nodeText(q.emitter)
            + ", model=" + q.model_name);
    }

    if (system.variable_names.empty()) {
        throw std::runtime_error("empty nonlinear MNA system");
    }

    return system;
}

Vector evaluateResidual(
    const Circuit& circuit, const NonlinearMnaSystem& system, const Vector& x)
{
    const auto node_index = makeNodeIndex(system.node_labels);
    const std::size_t node_count = system.node_labels.size();
    const std::size_t n = system.variable_names.size();
    if (x.size() != n) {
        throw std::runtime_error("nonlinear residual dimension mismatch");
    }

    Vector f(n, 0.0);

    for (const auto& r : circuit.resistors) {
        const double va = nodeVoltage(x, node_index, r.a);
        const double vb = nodeVoltage(x, node_index, r.b);
        const double current = (va - vb) / r.resistance;
        addResidualIfNode(f, node_index, r.a, current);
        addResidualIfNode(f, node_index, r.b, -current);
    }

    for (const auto& source : circuit.current_sources) {
        addResidualIfNode(f, node_index, source.p, source.current);
        addResidualIfNode(f, node_index, source.n, -source.current);
    }

    for (std::size_t k = 0; k < circuit.voltage_sources.size(); ++k) {
        const auto& source = circuit.voltage_sources[k];
        const std::size_t branch_col = node_count + k;
        const double branch_current = x.at(branch_col);
        addResidualIfNode(f, node_index, source.p, branch_current);
        addResidualIfNode(f, node_index, source.n, -branch_current);

        const double vp = nodeVoltage(x, node_index, source.p);
        const double vn = nodeVoltage(x, node_index, source.n);
        f.at(branch_col) += vp - vn - source.voltage;
    }

    for (const auto& q : circuit.bjts) {
        const auto& model = findBjtModel(circuit, q.model_name);
        const double vc = nodeVoltage(x, node_index, q.collector);
        const double vb = nodeVoltage(x, node_index, q.base);
        const double ve = nodeVoltage(x, node_index, q.emitter);
        const auto stamp = evaluateBjtStamp(model, vc, vb, ve);

        addResidualIfNode(f, node_index, q.collector, stamp.ic);
        addResidualIfNode(f, node_index, q.base, stamp.ib);
        addResidualIfNode(f, node_index, q.emitter, stamp.ie);
    }

    return f;
}

Matrix evaluateJacobian(
    const Circuit& circuit, const NonlinearMnaSystem& system, const Vector& x)
{
    const auto node_index = makeNodeIndex(system.node_labels);
    const std::size_t node_count = system.node_labels.size();
    const std::size_t n = system.variable_names.size();
    if (x.size() != n) {
        throw std::runtime_error("nonlinear Jacobian dimension mismatch");
    }

    Matrix jacobian(n, Vector(n, 0.0));

    for (const auto& r : circuit.resistors) {
        const double g = 1.0 / r.resistance;
        addJacobianNodeColumnIfNode(jacobian, node_index, r.a, r.a, g);
        addJacobianNodeColumnIfNode(jacobian, node_index, r.a, r.b, -g);
        addJacobianNodeColumnIfNode(jacobian, node_index, r.b, r.a, -g);
        addJacobianNodeColumnIfNode(jacobian, node_index, r.b, r.b, g);
    }

    for (std::size_t k = 0; k < circuit.voltage_sources.size(); ++k) {
        const auto& source = circuit.voltage_sources[k];
        const std::size_t branch_col = node_count + k;
        addJacobianExplicitColumnIfNode(
            jacobian, node_index, source.p, branch_col, 1.0);
        addJacobianExplicitColumnIfNode(
            jacobian, node_index, source.n, branch_col, -1.0);

        if (source.p != 0) {
            jacobian.at(branch_col).at(node_index.at(source.p)) += 1.0;
        }
        if (source.n != 0) {
            jacobian.at(branch_col).at(node_index.at(source.n)) -= 1.0;
        }
    }

    for (const auto& q : circuit.bjts) {
        const auto& model = findBjtModel(circuit, q.model_name);
        const double vc = nodeVoltage(x, node_index, q.collector);
        const double vb = nodeVoltage(x, node_index, q.base);
        const double ve = nodeVoltage(x, node_index, q.emitter);
        const auto stamp = evaluateBjtStamp(model, vc, vb, ve);

        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.collector, q.collector, stamp.dic_dvc);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.collector, q.base, stamp.dic_dvb);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.collector, q.emitter, stamp.dic_dve);

        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.base, q.collector, stamp.dib_dvc);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.base, q.base, stamp.dib_dvb);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.base, q.emitter, stamp.dib_dve);

        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.emitter, q.collector, stamp.die_dvc);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.emitter, q.base, stamp.die_dvb);
        addJacobianNodeColumnIfNode(
            jacobian, node_index, q.emitter, q.emitter, stamp.die_dve);
    }

    return jacobian;
}

void printNonlinearSystem(std::ostream& os, const NonlinearMnaSystem& system)
{
    os << "Variable order x:\n";
    for (std::size_t i = 0; i < system.variable_names.size(); ++i) {
        os << "  x[" << i << "] = " << system.variable_names[i] << '\n';
    }

    os << "\nNonlinear stamp log:\n";
    for (const auto& line : system.stamp_log) {
        os << "  - " << line << '\n';
    }
}

void printNamedVector(
    std::ostream& os,
    const std::vector<std::string>& names,
    const Vector& values,
    const std::string& title)
{
    os << '\n' << title << ":\n";
    os << std::scientific << std::setprecision(9);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const std::string name = i < names.size() ? names[i] : "x[" + std::to_string(i) + "]";
        os << "  " << std::setw(10) << name << " = " << std::setw(16)
           << values[i] << '\n';
    }
}

} // namespace dcsolve
