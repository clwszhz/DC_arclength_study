#include "dcsolve/Circuit.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dcsolve {
namespace {

std::string trim(const std::string& text)
{
    std::size_t first = 0;
    while (first < text.size()
           && std::isspace(static_cast<unsigned char>(text[first]))) {
        ++first;
    }

    std::size_t last = text.size();
    while (last > first
           && std::isspace(static_cast<unsigned char>(text[last - 1]))) {
        --last;
    }

    return text.substr(first, last - first);
}

std::string upperCopy(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return text;
}

std::vector<std::string> splitWords(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> words;
    std::string word;

    while (stream >> word) {
        words.push_back(word);
    }

    return words;
}

int parseNode(const std::string& text, int line_number)
{
    std::size_t consumed = 0;
    int value = 0;

    try {
        value = std::stoi(text, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": node must be an integer: " + text);
    }

    if (consumed != text.size() || value < 0) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": node must be a non-negative integer: "
                                 + text);
    }

    return value;
}

double parseNumber(const std::string& text, int line_number)
{
    std::size_t consumed = 0;
    double value = 0.0;

    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": value must be a plain number: " + text);
    }

    if (consumed != text.size()) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": suffixes like 1k/2meg are not supported yet: "
                                 + text);
    }

    return value;
}

void requireWordCount(
    const std::vector<std::string>& words, std::size_t expected, int line_number)
{
    if (words.size() != expected) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": expected " + std::to_string(expected)
                                 + " fields, got "
                                 + std::to_string(words.size()));
    }
}

std::string cleanParameterToken(std::string token)
{
    while (!token.empty() && token.front() == '(') {
        token.erase(token.begin());
    }
    while (!token.empty() && (token.back() == ')' || token.back() == ',')) {
        token.pop_back();
    }
    return token;
}

std::pair<std::string, std::string> splitAssignment(
    const std::string& token, int line_number)
{
    const auto pos = token.find('=');
    if (pos == std::string::npos || pos == 0 || pos + 1 == token.size()) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": model parameter must look like KEY=value: "
                                 + token);
    }

    return {upperCopy(token.substr(0, pos)), token.substr(pos + 1)};
}

bool hasBjtModel(const Circuit& circuit, const std::string& model_name)
{
    for (const auto& model : circuit.bjt_models) {
        if (model.name == model_name) {
            return true;
        }
    }
    return false;
}

BjtModel parseBjtModel(const std::vector<std::string>& words, int line_number)
{
    if (words.size() < 3) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": .MODEL needs at least name and type");
    }

    BjtModel model;
    model.name = upperCopy(words[1]);

    const std::string type = upperCopy(words[2]);
    if (type != "NPN") {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": only NPN BJT models are supported");
    }
    model.type = BjtType::Npn;

    for (std::size_t i = 3; i < words.size(); ++i) {
        const std::string token = cleanParameterToken(words[i]);
        if (token.empty()) {
            continue;
        }

        const auto [key, raw_value] = splitAssignment(token, line_number);
        const double value = parseNumber(raw_value, line_number);

        if (key == "IS") {
            model.saturation_current = value;
        } else if (key == "AF") {
            model.alpha_forward = value;
        } else if (key == "AR") {
            model.alpha_reverse = value;
        } else if (key == "BF") {
            model.alpha_forward = value / (value + 1.0);
        } else if (key == "BR") {
            model.alpha_reverse = value / (value + 1.0);
        } else if (key == "N") {
            model.exponential_coefficient = value;
        } else {
            throw std::runtime_error("line " + std::to_string(line_number)
                                     + ": unsupported BJT model parameter: "
                                     + key);
        }
    }

    if (model.saturation_current <= 0.0) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": BJT IS must be positive");
    }
    if (model.alpha_forward <= 0.0 || model.alpha_forward > 1.0
        || model.alpha_reverse <= 0.0 || model.alpha_reverse > 1.0) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": BJT AF/AR must be in (0, 1]");
    }
    if (model.exponential_coefficient <= 0.0) {
        throw std::runtime_error("line " + std::to_string(line_number)
                                 + ": BJT N must be positive");
    }

    return model;
}

} // namespace

Circuit parseNetlistFile(const std::string& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open netlist: " + path);
    }

    Circuit circuit;
    std::string line;
    int line_number = 0;

    while (std::getline(input, line)) {
        ++line_number;

        const auto comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        line = trim(line);
        if (line.empty() || line[0] == '*') {
            continue;
        }

        const auto words = splitWords(line);
        const std::string keyword = upperCopy(words[0]);

        if (keyword == ".OP") {
            requireWordCount(words, 1, line_number);
            circuit.has_op = true;
            continue;
        }

        if (keyword == ".MODEL") {
            const auto model = parseBjtModel(words, line_number);
            if (hasBjtModel(circuit, model.name)) {
                throw std::runtime_error("line " + std::to_string(line_number)
                                         + ": duplicate model name: "
                                         + model.name);
            }
            circuit.bjt_models.push_back(model);
            continue;
        }

        if (keyword == ".START") {
            if (words.size() < 2) {
                throw std::runtime_error("line " + std::to_string(line_number)
                                         + ": .START needs at least one value");
            }
            circuit.initial_guess.clear();
            for (std::size_t i = 1; i < words.size(); ++i) {
                circuit.initial_guess.push_back(parseNumber(words[i], line_number));
            }
            continue;
        }

        if (keyword == ".END") {
            break;
        }

        const char type = static_cast<char>(
            std::toupper(static_cast<unsigned char>(words[0][0])));

        switch (type) {
        case 'R': {
            requireWordCount(words, 4, line_number);
            Resistor r;
            r.name = words[0];
            r.a = parseNode(words[1], line_number);
            r.b = parseNode(words[2], line_number);
            r.resistance = parseNumber(words[3], line_number);
            if (r.resistance <= 0.0) {
                throw std::runtime_error("line " + std::to_string(line_number)
                                         + ": resistance must be positive");
            }
            circuit.resistors.push_back(r);
            break;
        }

        case 'I': {
            requireWordCount(words, 4, line_number);
            CurrentSource source;
            source.name = words[0];
            source.p = parseNode(words[1], line_number);
            source.n = parseNode(words[2], line_number);
            source.current = parseNumber(words[3], line_number);
            circuit.current_sources.push_back(source);
            break;
        }

        case 'V': {
            requireWordCount(words, 4, line_number);
            VoltageSource source;
            source.name = words[0];
            source.p = parseNode(words[1], line_number);
            source.n = parseNode(words[2], line_number);
            source.voltage = parseNumber(words[3], line_number);
            circuit.voltage_sources.push_back(source);
            break;
        }

        case 'Q': {
            requireWordCount(words, 5, line_number);
            Bjt q;
            q.name = words[0];
            q.collector = parseNode(words[1], line_number);
            q.base = parseNode(words[2], line_number);
            q.emitter = parseNode(words[3], line_number);
            q.model_name = upperCopy(words[4]);
            circuit.bjts.push_back(q);
            break;
        }

        default:
            throw std::runtime_error("line " + std::to_string(line_number)
                                     + ": unsupported element or command: "
                                     + words[0]);
        }
    }

    if (!circuit.has_op) {
        throw std::runtime_error("netlist must contain .OP for this first demo");
    }

    for (const auto& q : circuit.bjts) {
        if (!hasBjtModel(circuit, q.model_name)) {
            throw std::runtime_error("BJT " + q.name
                                     + " references undefined model: "
                                     + q.model_name);
        }
    }

    return circuit;
}

} // namespace dcsolve
