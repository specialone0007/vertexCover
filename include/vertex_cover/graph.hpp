#pragma once
// Undirected simple graph stored as adjacency lists. Vertices are 0..V-1.
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

namespace vc {

using Edge = std::pair<int, int>;

class Graph {
public:
    explicit Graph(int vertexCount);

    int vertexCount() const { return static_cast<int>(adj_.size()); }
    int edgeCount() const { return edgeCount_; }
    const std::vector<int>& neighbors(int v) const { return adj_[v]; }

    // Adds an undirected edge. Self-loops and duplicates are rejected (returns false).
    bool addEdge(int u, int v);
    bool hasEdge(int u, int v) const;

    // Every edge exactly once, with first < second.
    std::vector<Edge> edges() const;

    // Random simple graph with exactly `edgeCount` edges where every vertex has
    // degree >= 1 (when edgeCount permits). Throws std::invalid_argument if
    // edgeCount exceeds V*(V-1)/2.
    static Graph random(int vertexCount, int edgeCount, std::mt19937& rng);

private:
    std::vector<std::vector<int>> adj_;
    int edgeCount_ = 0;
};

}  // namespace vc
