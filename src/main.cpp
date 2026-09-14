// vertex_cover CLI
//   solve   <V> <E> [--seed S] [--exact]         greedy (and optionally exact) cover of one random graph
//   solve   --file <edges.txt> [--exact]         same, on a graph read from an edge list ("u v" per line, 0-based)
//   verify  <graphs> <Vmin> <Vmax> <Emax>        check greedy output is a valid cover on random graphs
//   quality <trials>                             greedy vs exact ratio on small graphs (V = 5..20)
//   bench   <iterations> [--out DIR]             O(V+E) running-time study, one CSV per configuration
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "vertex_cover/algorithms.hpp"
#include "vertex_cover/graph.hpp"
#include "vertex_cover/stats.hpp"

namespace {

int usage() {
    std::cerr <<
        "usage:\n"
        "  vertex_cover solve   <V> <E> [--seed S] [--exact]\n"
        "  vertex_cover solve   --file <edges.txt> [--exact]\n"
        "  vertex_cover verify  <graphs> <Vmin> <Vmax> <Emax> [--seed S]\n"
        "  vertex_cover quality <trials> [--seed S]\n"
        "  vertex_cover bench   <iterations> [--out DIR] [--seed S]\n";
    return 2;
}

struct Args {
    std::vector<std::string> positional;
    unsigned seed = 42;
    bool exact = false;
    std::string out = ".";
    std::string file;
};

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--seed" && i + 1 < argc) a.seed = static_cast<unsigned>(std::stoul(argv[++i]));
        else if (s == "--out" && i + 1 < argc) a.out = argv[++i];
        else if (s == "--file" && i + 1 < argc) a.file = argv[++i];
        else if (s == "--exact") a.exact = true;
        else a.positional.push_back(s);
    }
    return a;
}

double timeGreedyMs(const vc::Graph& g) {
    // One untimed warm-up run so page faults and cache misses from building the
    // graph do not land inside the measured window.
    volatile std::size_t warm = vc::greedyCover(g).size();
    (void)warm;
    const auto t0 = std::chrono::steady_clock::now();
    volatile std::size_t sink = vc::greedyCover(g).size();
    (void)sink;
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void writeRow(std::ostream& os, int size, const vc::Summary& s) {
    os << size << ',' << s.mean << ',' << s.stddev << ',' << s.stderror << ','
       << s.ci90Low << ',' << s.ci90High << ',' << s.ci95Low << ',' << s.ci95High << '\n';
}

// Edge list: one "u v" pair per line, 0-based; '#' starts a comment. V = max vertex id + 1.
vc::Graph readEdgeList(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<vc::Edge> edges;
    int maxV = -1;
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ls(line);
        int u, v;
        if (!(ls >> u >> v)) continue;
        if (u < 0 || v < 0) throw std::runtime_error("negative vertex id in " + path);
        edges.emplace_back(u, v);
        maxV = std::max({maxV, u, v});
    }
    vc::Graph g(maxV + 1);
    for (auto [u, v] : edges) g.addEdge(u, v);
    return g;
}

int cmdSolve(const Args& a) {
    vc::Graph g(0);
    if (!a.file.empty()) {
        if (!a.positional.empty()) return usage();
        g = readEdgeList(a.file);
    } else {
        if (a.positional.size() != 2) return usage();
        std::mt19937 rng(a.seed);
        g = vc::Graph::random(std::stoi(a.positional[0]), std::stoi(a.positional[1]), rng);
    }
    const auto greedy = vc::greedyCover(g);
    std::cout << "graph: V=" << g.vertexCount() << " E=" << g.edgeCount() << "\n";
    std::cout << "greedy cover (" << greedy.size() << "):";
    for (int v : greedy) std::cout << ' ' << v;
    std::cout << "\nvalid: " << (vc::isVertexCover(g, greedy) ? "yes" : "no") << "\n";
    if (a.exact) {
        const auto opt = vc::exactCover(g);
        std::cout << "exact cover (" << opt.size() << "):";
        for (int v : opt) std::cout << ' ' << v;
        std::cout << "\nratio: " << static_cast<double>(greedy.size()) / static_cast<double>(opt.size()) << "\n";
    }
    return 0;
}

int cmdVerify(const Args& a) {
    if (a.positional.size() != 4) return usage();
    const int graphs = std::stoi(a.positional[0]);
    const int vMin = std::stoi(a.positional[1]), vMax = std::stoi(a.positional[2]);
    const int eMax = std::stoi(a.positional[3]);
    std::mt19937 rng(a.seed);
    std::uniform_int_distribution<int> pickV(vMin, vMax);
    int ok = 0;
    for (int i = 0; i < graphs; ++i) {
        const int V = pickV(rng);
        const long long cap = std::min<long long>(eMax, 1LL * V * (V - 1) / 2);
        std::uniform_int_distribution<long long> pickE(std::min<long long>(V, cap), cap);
        const auto g = vc::Graph::random(V, static_cast<int>(pickE(rng)), rng);
        if (vc::isVertexCover(g, vc::greedyCover(g))) ++ok;
    }
    std::cout << "valid covers: " << ok << "/" << graphs << "\n";
    return ok == graphs ? 0 : 1;
}

int cmdQuality(const Args& a) {
    if (a.positional.size() != 1) return usage();
    const int trials = std::stoi(a.positional[0]);
    std::mt19937 rng(a.seed);
    std::cout << "V,trials,mean_quality(opt/greedy),worst_ratio(greedy/opt),optimal_hits\n";
    for (int V : {5, 6, 7, 8, 9, 10, 12, 15, 18, 20}) {
        std::uniform_int_distribution<int> pickE(V - 1, V * (V - 1) / 2);
        double sumQ = 0, worst = 1;
        int hits = 0;
        for (int t = 0; t < trials; ++t) {
            const auto g = vc::Graph::random(V, pickE(rng), rng);
            const double h = static_cast<double>(vc::greedyCover(g).size());
            const double o = static_cast<double>(vc::exactCover(g).size());
            sumQ += o / h;
            worst = std::max(worst, h / o);
            if (h == o) ++hits;
        }
        std::cout << V << ',' << trials << ',' << sumQ / trials << ',' << worst << ',' << hits << '\n';
    }
    return 0;
}

int cmdBench(const Args& a) {
    if (a.positional.size() != 1) return usage();
    const int iters = std::stoi(a.positional[0]);
    std::mt19937 rng(a.seed);
    const char* header = "size,mean_ms,stddev,stderr,ci90_low,ci90_high,ci95_low,ci95_high\n";

    // Sizes are large enough that one greedy run takes milliseconds, so the
    // measurement is not dominated by clock resolution.
    {   // E fixed at 20,000; V = 20k .. 200k
        std::ofstream f(a.out + "/edges-fixed-" + std::to_string(iters) + "-iter.csv");
        f << header;
        for (int V = 20000; V <= 200000; V += 20000) {
            std::vector<double> t;
            for (int i = 0; i < iters; ++i) t.push_back(timeGreedyMs(vc::Graph::random(V, 20000, rng)));
            writeRow(f, V, vc::summarize(t));
        }
    }
    {   // V fixed at 20,000; E = 20k .. 1M
        std::ofstream f(a.out + "/vertices-fixed-" + std::to_string(iters) + "-iter.csv");
        f << header;
        for (int E = 20000; E <= 1000000; E += 140000) {
            std::vector<double> t;
            for (int i = 0; i < iters; ++i) t.push_back(timeGreedyMs(vc::Graph::random(20000, E, rng)));
            writeRow(f, E, vc::summarize(t));
        }
    }
    std::cout << "wrote 2 CSV files to " << a.out << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return usage();
    const std::string cmd = argv[1];
    const Args a = parse(argc, argv);
    try {
        if (cmd == "solve") return cmdSolve(a);
        if (cmd == "verify") return cmdVerify(a);
        if (cmd == "quality") return cmdQuality(a);
        if (cmd == "bench") return cmdBench(a);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return usage();
}
