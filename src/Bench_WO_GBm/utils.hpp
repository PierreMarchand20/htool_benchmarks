#ifndef UTILS_HPP
#define UTILS_HPP
#include <cmath>
#include <numeric>

void compute_standard_deviation(const double arr[], int size, double &mean, double &standard_deviation) {
    double sum                 = std::accumulate(arr, arr + size, 0.0);
    mean                       = sum / size;
    double squared_differences = std::inner_product(arr, arr + size, arr, 0.0, std::plus<double>(), [mean](double x, double y) { return (x - mean) * (x - mean); });
    double variance            = squared_differences / (size - 1);
    standard_deviation         = std::sqrt(variance);
}

std::vector<int> geometric_progression(int a, int r, int n) {
    // Computes the sequence of n terms of a geometric progression.
    //
    // Parameters:
    // a: The first term of the progression.
    // r: The common ratio of the progression.
    // n: The number of terms to compute.
    //
    // Returns:
    // A vector containing the sequence of n terms of the progression.
    std::vector<int> sequence(n);
    sequence[0] = a;
    for (int i = 1; i < n; i++) {
        sequence[i] = sequence[i - 1] * r;
    }
    return sequence;
}

#endif