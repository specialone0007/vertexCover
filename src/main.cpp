// vertex_cover CLI
//   solve   <V> <E> [--seed S] [--exact]         greedy (and optionally exact) cover of one random graph
//   verify  <graphs> <Vmin> <Vmax> <Emax>        check greedy output is a valid cover on random graphs
//   quality <trials>                             greedy vs exact ratio on small graphs (V = 5..20)
//   bench   <iterations> [--out DIR]             O(V+E) running-time study, one CSV per configuration
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
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
};

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--seed" && i + 1 < argc) a.seed = static_cast<unsigned>(std::stoul(argv[++i]));
        else if (s == "--out" && i + 1 < argc) a.out = argv[++i];
        else if (s == "--exact") a.exact = true;
        else a.positional.push_back(s);
    }
    return a;
}

double timeGreedyMs(const vc::Graph& g) {
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

int cmdSolve(const Args& a) {
    if (a.positional.size() != 2) return usage();
    const int V = std::stoi(a.positional[0]), E = std::stoi(a.positional[1]);
    std::mt19937 rng(a.seed);
    const auto g = vc::Graph::random(V, E, rng);
    const auto greedy = vc::greedyCover(g);
    std::cout << "graph: V=" << V << " E=" << g.edgeCount() << "\n";
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

    {   // E fixed at 200, V = 100..1000
        std::ofstream f(a.out + "/edges-fixed-200-" + std::to_string(iters) + "-iter.csv");
        f << header;
        for (int V = 100; V <= 1000; V += 100) {
            std::vector<double> t;
            for (int i = 0; i < iters; ++i) t.push_back(timeGreedyMs(vc::Graph::random(V, 200, rng)));
            writeRow(f, V, vc::summarize(t));
        }
    }
    {   // V fixed at 200, E = 200..4700
        std::ofstream f(a.out + "/vertices-fixed-200-" + std::to_string(iters) + "-iter.csv");
        f << header;
        for (int E = 200; E <= 4700; E += 500) {
            std::vector<double> t;
            for (int i = 0; i < iters; ++i) t.push_back(timeGreedyMs(vc::Graph::random(200, E, rng)));
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
