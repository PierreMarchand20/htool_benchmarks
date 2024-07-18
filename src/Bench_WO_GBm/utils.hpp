#include <cmath>
#include <numeric>

void compute_standard_deviation(const double arr[], int size, double &mean, double &standard_deviation) {
    double sum                 = std::accumulate(arr, arr + size, 0.0);
    mean                       = sum / size;
    double squared_differences = std::inner_product(arr, arr + size, arr, 0.0, std::plus<double>(), [mean](double x, double y) { return (x - mean) * (x - mean); });
    double variance            = squared_differences / (size - 1);
    standard_deviation         = std::sqrt(variance);
}