# vertex-cover

[![ci](https://github.com/specialone0007/vertexCover/actions/workflows/ci.yml/badge.svg)](https://github.com/specialone0007/vertexCover/actions/workflows/ci.yml)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)
![license](https://img.shields.io/badge/license-MIT-green)

A greedy **2-approximation for Minimum Vertex Cover**, an exact brute-force solver to
measure it against, and a small experimental harness: running-time study with confidence
intervals, correctness fuzzing, and approximation-quality sampling.

Started as a Sabancı University CS301 (Algorithms) term project in 2019, rewritten in 2026
as a proper library + CLI with tests and CI.

**[Read the report (PDF, 8 pages)](docs/report.pdf)** — NP-completeness proof by reduction from
3-SAT, the 2-approximation proof, fresh experiments, and a postmortem of the bugs in the 2019 code.
Source in [`docs/report.md`](docs/report.md); the original 2019 group report is kept in
[`docs/legacy/`](docs/legacy/).

![Greedy vs optimum on a star and on the Petersen graph](docs/figures/greedy-vs-optimum.png)

## The problem

Given an undirected graph *G = (V, E)*, a **vertex cover** is a set *C ⊆ V* such that every
edge has at least one endpoint in *C*. Finding a **minimum** vertex cover is NP-complete, so
in practice we trade optimality for a polynomial-time guarantee.

## The algorithm

Maximal-matching greedy:

```
C ← ∅
while there is an edge (u, v) with u ∉ C and v ∉ C:
    C ← C ∪ {u, v}
return C
```

- **Time:** O(V + E). Each adjacency list is scanned at most once.
- **Guarantee:** |C| ≤ 2 · |OPT|. The chosen edges form a matching *A*, so |C| = 2|A|. Any
  cover, including the optimum, must contain at least one endpoint of every edge in *A*,
  and those endpoints are distinct, so |OPT| ≥ |A|. Hence |C| ≤ 2 |OPT|.
- **Tightness:** a star graph forces greedy to pick 2 vertices where 1 suffices. The bound is
  reached in practice as well (see *Results*).

The exact solver enumerates vertex subsets in order of increasing size as 64-bit masks
(Gosper's hack), stopping at the first subset that covers every edge. It exists purely as a
ground truth for graphs up to ~24 vertices.

## Layout

```
include/vertex_cover/   public headers: Graph, greedyCover, exactCover, isVertexCover, summarize
src/                    implementation + CLI (main.cpp)
tests/                  dependency-free unit tests (ctest)
examples/               edge-list graphs (Petersen, star) for `solve --file`
benchmarks/             results/ CSVs from `vertex_cover bench`, plot.py renders them
docs/                   report.md → report.pdf, figures/, legacy/ (2019 report)
```

## Build

Requires CMake ≥ 3.16 and any C++17 compiler (GCC, Clang, MSVC).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Usage

```bash
# one random graph, greedy vs exact
./build/vertex_cover solve 12 20 --exact
# graph: V=12 E=20
# greedy cover (10): 0 4 1 9 2 11 3 7 5 8
# valid: yes
# exact cover (7): 0 2 3 7 8 9 10
# ratio: 1.42857

# any graph from an edge list ("u v" per line, 0-based, # comments)
./build/vertex_cover solve --file examples/petersen.txt --exact
# graph: V=10 E=15
# greedy cover (10): 0 1 2 3 4 9 5 7 6 8
# exact cover (6): 1 3 4 5 6 7
# ratio: 1.66667

# fuzz: is greedy output always a valid cover?
./build/vertex_cover verify 10000 1000 10000 20000

# approximation quality on small graphs (needs exact solver, V <= 20)
./build/vertex_cover quality 500

# running-time study; writes two CSVs (E fixed / V fixed), then plot them
./build/vertex_cover bench 50 --out benchmarks/results
python benchmarks/plot.py
```

All commands take `--seed <n>` (default 42) so runs are reproducible.

## Results

**Correctness.** 10,000 random graphs with V ∈ [1000, 10000], E ∈ [V, 20000]: greedy returned
a valid cover every time. Also enforced by a ctest (`greedy_is_always_a_cover`).

**Quality** = |OPT| / |greedy| ∈ [0.5, 1], 500 random graphs per size, seed 42:

| V  | mean quality | worst ratio | exact hits / 500 |
|----|--------------|-------------|------------------|
| 5  | 0.74 | 2.0 | 81 |
| 8  | 0.70 | 2.0 | 16 |
| 10 | 0.72 | 2.0 | 14 |
| 15 | 0.77 | 2.0 | 9 |
| 20 | 0.78 | 2.0 | 2 |

![quality](docs/figures/quality.png)

Greedy typically lands 25–30 % above optimum, hits the theoretical worst case (ratio 2) at
every size, and almost never finds the exact optimum once V > 10.

**Running time.** 50 random graphs per size, warm-up run per graph, 95 % CI bands.

| E fixed at 20,000 · V = 20k…200k | V fixed at 20,000 · E = 20k…1M |
|---|---|
| ![](docs/figures/time-edges-fixed.png) | ![](docs/figures/time-vertices-fixed.png) |

Linear in V. Nearly flat in E: the inner loop stops at the first uncovered neighbour, and on
dense random graphs almost every vertex is covered early, so O(V + E) is a worst-case bound
that random inputs never approach. Discussion in the report, §5.3.

## What changed since 2019

The original was a single 400-line `main.cpp` with three `cin`-driven modes. The rewrite
keeps the algorithms and experiment design but fixes several real bugs found on the way:

- The timed greedy variant marked *all neighbours* of the chosen edge as visited, so it
  could skip edges and return a set that was **not a vertex cover**. It passed the 2019
  correctness check because…
- …`isVertexCover` tested *domination* (every vertex is in or adjacent to the set), not edge
  coverage. Both are now correct and cross-checked in tests.
- The exact solver mixed 1-based combinations with 0-based adjacency and rebuilt a
  neighbour list per subset. It is now a bitmask search with early exit.
- Random graph generation could spin forever when E exceeded V(V−1)/2 and could emit
  self-loops while patching isolated vertices. Both now rejected.
- `system("pause")`, raw `new[]` without `delete[]`, and `using namespace std` in headers
  are gone; everything is RAII and lives in namespace `vc`.

## License

MIT. Original project by Furkan Reha Tutaş, Tugay Garib, Derya Bensu Çakar, Meltem Arslan
and Emre Hilmi Songur (CS301, Fall 2019, instructor Hüsnü Yenigün). 2026 rewrite by
Furkan Reha Tutaş.
