#pragma once
// Descriptive statistics for benchmark samples.
#include <vector>

namespace vc {

struct Summary {
    double mean = 0;
    double stddev = 0;       // population standard deviation
    double stderror = 0;     // stddev / sqrt(n)
    double ci90Low = 0, ci90High = 0;  // z = 1.645
    double ci95Low = 0, ci95High = 0;  // z = 1.960
};

Summary summarize(const std::vector<double>& samples);

}  // namespace vc
