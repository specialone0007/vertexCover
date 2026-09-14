// Dependency-free unit tests. Non-zero exit code on failure (consumed by ctest).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>

#include "vertex_cover/algorithms.hpp"
#include "vertex_cover/graph.hpp"
#include "vertex_cover/stats.hpp"

static int failures = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            ++failures;                                                                \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " << #cond << "\n"; \
        }                                                                              \
    } while (0)

static void testGraphBasics() {
    vc::Graph g(4);
    CHECK(g.addEdge(0, 1));
    CHECK(!g.addEdge(1, 0));   // duplicate
    CHECK(!g.addEdge(2, 2));   // self-loop
    CHECK(g.addEdge(2, 3));
    CHECK(g.edgeCount() == 2);
    CHECK(g.edges().size() == 2);
    CHECK(g.hasEdge(3, 2));
}

static void testStar() {  // center alone is optimal; greedy takes one edge = 2 vertices
    vc::Graph g(6);
    for (int v = 1; v < 6; ++v) g.addEdge(0, v);
    const auto greedy = vc::greedyCover(g);
    const auto exact = vc::exactCover(g);
    CHECK(vc::isVertexCover(g, greedy));
    CHECK(greedy.size() == 2);
    CHECK(exact.size() == 1 && exact[0] == 0);
}

static void testPath() {  // 0-1-2-3: optimum {1,2}
    vc::Graph g(4);
    g.addEdge(0, 1); g.addEdge(1, 2); g.addEdge(2, 3);
    CHECK(vc::exactCover(g).size() == 2);
    CHECK(vc::isVertexCover(g, vc::greedyCover(g)));
    CHECK(!vc::isVertexCover(g, {0, 3}));
    CHECK(vc::isVertexCover(g, {1, 2}));
}

static void testEmptyAndEdgeless() {
    vc::Graph g0(0);
    CHECK(vc::greedyCover(g0).empty());
    CHECK(vc::exactCover(g0).empty());
    vc::Graph g3(3);
    CHECK(vc::greedyCover(g3).empty());
    CHECK(vc::isVertexCover(g3, {}));
}

static void testRandomInvariants() {
    std::mt19937 rng(7);
    for (int t = 0; t < 300; ++t) {
        std::uniform_int_distribution<int> pickV(2, 14);
        const int V = pickV(rng);
        std::uniform_int_distribution<int> pickE(1, V * (V - 1) / 2);
        const auto g = vc::Graph::random(V, pickE(rng), rng);
        const auto greedy = vc::greedyCover(g);
        const auto exact = vc::exactCover(g);
        CHECK(vc::isVertexCover(g, greedy));
        CHECK(vc::isVertexCover(g, exact));
        CHECK(exact.size() <= greedy.size());
        CHECK(greedy.size() <= 2 * exact.size());   // 2-approximation bound
        CHECK(greedy.size() % 2 == 0);               // matching endpoints come in pairs
    }
}

static void testRandomGraphShape() {
    std::mt19937 rng(1);
    const auto g = vc::Graph::random(50, 120, rng);
    CHECK(g.edgeCount() == 120);
    for (int v = 0; v < 50; ++v) CHECK(!g.neighbors(v).empty());
    bool threw = false;
    try { vc::Graph::random(4, 7, rng); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

static void testStats() {
    const auto s = vc::summarize({2, 4, 4, 4, 5, 5, 7, 9});
    CHECK(std::abs(s.mean - 5.0) < 1e-12);
    CHECK(std::abs(s.stddev - 2.0) < 1e-12);
    CHECK(s.ci95Low < s.mean && s.mean < s.ci95High);
}

int main() {
    testGraphBasics();
    testStar();
    testPath();
    testEmptyAndEdgeless();
    testRandomInvariants();
    testRandomGraphShape();
    testStats();
    if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
    std::cout << "all tests passed\n";
    return 0;
}
