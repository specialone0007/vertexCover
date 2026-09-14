#include "vertex_cover/graph.hpp"

#include <algorithm>
#include <stdexcept>

namespace vc {

Graph::Graph(int vertexCount) : adj_(vertexCount) {
    if (vertexCount < 0) throw std::invalid_argument("vertexCount must be >= 0");
}

bool Graph::hasEdge(int u, int v) const {
    const auto& list = adj_[u];
    return std::find(list.begin(), list.end(), v) != list.end();
}

bool Graph::addEdge(int u, int v) {
    if (u == v || hasEdge(u, v)) return false;
    adj_[u].push_back(v);
    adj_[v].push_back(u);
    ++edgeCount_;
    return true;
}

std::vector<Edge> Graph::edges() const {
    std::vector<Edge> out;
    out.reserve(edgeCount_);
    for (int u = 0; u < vertexCount(); ++u)
        for (int v : adj_[u])
            if (u < v) out.emplace_back(u, v);
    return out;
}

Graph Graph::random(int vertexCount, int edgeCount, std::mt19937& rng) {
    const long long maxEdges = 1LL * vertexCount * (vertexCount - 1) / 2;
    if (edgeCount < 0 || edgeCount > maxEdges)
        throw std::invalid_argument("edgeCount out of range for a simple graph");

    Graph g(vertexCount);
    if (vertexCount < 2) return g;
    std::uniform_int_distribution<int> pick(0, vertexCount - 1);

    // Guarantee degree >= 1 first when the budget allows: attach every isolated
    // vertex to a random other vertex.
    if (edgeCount >= (vertexCount + 1) / 2) {
        for (int u = 0; u < vertexCount && g.edgeCount() < edgeCount; ++u) {
            if (!g.adj_[u].empty()) continue;
            int v;
            do { v = pick(rng); } while (v == u);
            g.addEdge(u, v);
        }
    }
    while (g.edgeCount() < edgeCount) {
        g.addEdge(pick(rng), pick(rng));  // rejected if self-loop or duplicate
    }
    return g;
}

}  // namespace vc
