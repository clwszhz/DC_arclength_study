#pragma once

#include <string>
#include <vector>

namespace dcsolve {

enum class BjtType {
    Npn,
};

struct Resistor {
    std::string name;
    int a = 0;
    int b = 0;
    double resistance = 0.0;
};

struct CurrentSource {
    std::string name;
    int p = 0;
    int n = 0;
    double current = 0.0;
};

struct VoltageSource {
    std::string name;
    int p = 0;
    int n = 0;
    double voltage = 0.0;
};

struct BjtModel {
    std::string name;
    BjtType type = BjtType::Npn;
    double saturation_current = 1e-16;
    double alpha_forward = 0.99;
    double alpha_reverse = 0.5;
    double exponential_coefficient = 38.78;
};

struct Bjt {
    std::string name;
    int collector = 0;
    int base = 0;
    int emitter = 0;
    std::string model_name;
};

struct Circuit {
    std::vector<Resistor> resistors;
    std::vector<CurrentSource> current_sources;
    std::vector<VoltageSource> voltage_sources;
    std::vector<BjtModel> bjt_models;
    std::vector<Bjt> bjts;
    std::vector<double> initial_guess;
    bool has_op = false;

    std::vector<int> nodeLabels() const;
};

Circuit parseNetlistFile(const std::string& path);

} // namespace dcsolve
