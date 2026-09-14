#pragma once
// Vertex-cover solvers and checkers.
#include <vector>

#include "vertex_cover/graph.hpp"

namespace vc {

// Greedy 2-approximation via maximal matching:
// pick any uncovered edge (u, v), add both endpoints, repeat.
// Runs in O(V + E). |result| <= 2 * |optimum|.
std::vector<int> greedyCover(const Graph& g);

// Exact minimum vertex cover by enumerating vertex subsets in order of
// increasing size (bitmask). Exponential; intended for V <= ~24.
std::vector<int> exactCover(const Graph& g);

// True iff every edge has at least one endpoint in `cover`.
bool isVertexCover(const Graph& g, const std::vector<int>& cover);

}  // namespace vc
