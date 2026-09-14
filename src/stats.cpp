#include "vertex_cover/stats.hpp"

#include <cmath>

namespace vc {

Summary summarize(const std::vector<double>& s) {
    Summary out;
    if (s.empty()) return out;
    const double n = static_cast<double>(s.size());
    double sum = 0;
    for (double x : s) sum += x;
    out.mean = sum / n;
    double sq = 0;
    for (double x : s) sq += (x - out.mean) * (x - out.mean);
    out.stddev = std::sqrt(sq / n);
    out.stderror = out.stddev / std::sqrt(n);
    out.ci90Low = out.mean - 1.645 * out.stderror;
    out.ci90High = out.mean + 1.645 * out.stderror;
    out.ci95Low = out.mean - 1.960 * out.stderror;
    out.ci95High = out.mean + 1.960 * out.stderror;
    return out;
}

}  // namespace vc
