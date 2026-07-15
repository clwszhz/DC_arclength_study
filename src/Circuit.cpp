#include "dcsolve/Circuit.h"

#include <algorithm>
#include <set>

namespace dcsolve {

std::vector<int> Circuit::nodeLabels() const
{
    std::set<int> labels;

    auto add_node = [&labels](int node) {
        if (node != 0) {
            labels.insert(node);
        }
    };

    for (const auto& r : resistors) {
        add_node(r.a);
        add_node(r.b);
    }

    for (const auto& i : current_sources) {
        add_node(i.p);
        add_node(i.n);
    }

    for (const auto& v : voltage_sources) {
        add_node(v.p);
        add_node(v.n);
    }

    for (const auto& q : bjts) {
        add_node(q.collector);
        add_node(q.base);
        add_node(q.emitter);
    }

    return std::vector<int>(labels.begin(), labels.end());
}

} // namespace dcsolve
