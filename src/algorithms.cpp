#include "vertex_cover/algorithms.hpp"

#include <cstdint>
#include <stdexcept>

namespace vc {

std::vector<int> greedyCover(const Graph& g) {
    const int n = g.vertexCount();
    std::vector<char> inCover(n, 0);
    std::vector<int> cover;
    for (int u = 0; u < n; ++u) {
        if (inCover[u]) continue;
        for (int v : g.neighbors(u)) {
            if (inCover[v]) continue;
            // (u, v) is uncovered: take both endpoints, then stop scanning u.
            inCover[u] = inCover[v] = 1;
            cover.push_back(u);
            cover.push_back(v);
            break;
        }
    }
    return cover;
}

bool isVertexCover(const Graph& g, const std::vector<int>& cover) {
    std::vector<char> inCover(g.vertexCount(), 0);
    for (int v : cover) {
        if (v < 0 || v >= g.vertexCount()) return false;
        inCover[v] = 1;
    }
    for (int u = 0; u < g.vertexCount(); ++u) {
        if (inCover[u]) continue;
        for (int v : g.neighbors(u))
            if (!inCover[v]) return false;
    }
    return true;
}

namespace {
// Gosper's hack: next bitmask with the same popcount.
std::uint64_t nextSameWeight(std::uint64_t x) {
    const std::uint64_t c = x & -x;
    const std::uint64_t r = x + c;
    return (((r ^ x) >> 2) / c) | r;
}
}  // namespace

std::vector<int> exactCover(const Graph& g) {
    const int n = g.vertexCount();
    if (n > 62) throw std::invalid_argument("exactCover supports at most 62 vertices");
    const auto edges = g.edges();
    if (edges.empty()) return {};

    std::vector<std::uint64_t> edgeMask;
    edgeMask.reserve(edges.size());
    for (auto [u, v] : edges) edgeMask.push_back((1ULL << u) | (1ULL << v));

    const std::uint64_t limit = (n == 62) ? ~0ULL : ((1ULL << n) - 1);
    for (int k = 1; k <= n; ++k) {
        std::uint64_t mask = (1ULL << k) - 1;
        while (mask <= limit) {
            bool ok = true;
            for (std::uint64_t em : edgeMask)
                if ((em & mask) == 0) { ok = false; break; }
            if (ok) {
                std::vector<int> cover;
                for (int v = 0; v < n; ++v)
                    if (mask >> v & 1) cover.push_back(v);
                return cover;
            }
            const std::uint64_t next = nextSameWeight(mask);
            if (next <= mask) break;  // overflow
            mask = next;
        }
    }
    return {};  // unreachable: the full vertex set always covers
}

}  // namespace vc
