#ifndef UTILS_HPP
#define UTILS_HPP
#include <cmath>
#include <numeric>

/**
 * @brief Compute the mean and standard deviation of an array of double values.
 *
 * @param arr The array of double values.
 * @param size The number of elements in the array.
 * @param mean The mean of the array.
 * @param standard_deviation The standard deviation of the array.
 */
inline void compute_standard_deviation(const double arr[], int size, double &mean, double &standard_deviation) {
    // Compute the sum of all elements in the array.
    double sum = std::accumulate(arr, arr + size, 0.0);

    // Compute the mean of the array.
    mean = sum / size;

    // Compute the squared differences between each element and the mean.
    double squared_differences = std::inner_product(arr, arr + size, arr, 0.0, std::plus<double>(), [mean](double x, double y) { return (x - mean) * (x - mean); });

    // Compute the variance as the average of the squared differences.
    double variance = squared_differences / (size - 1);

    // Compute the standard deviation as the square root of the variance.
    standard_deviation = std::sqrt(variance);
}

/**
 * @brief Compute a geometric progression of integers.
 *
 * @param firstTerm The first term of the geometric progression.
 * @param commonRatio The common ratio of the geometric progression.
 * @param numTerms The number of terms in the geometric progression.
 * @return A vector of the first <code>numTerms</code> terms of the geometric progression.
 */
inline std::vector<int> geometricProgression(int firstTerm, int commonRatio, int numTerms) {
    std::vector<int> sequence(numTerms);
    // The first term of the sequence is the first term of the geometric progression.
    sequence[0] = firstTerm;
    // Compute the rest of the sequence by multiplying each term by the common ratio.
    for (int i = 1; i < numTerms; i++) {
        sequence[i] = sequence[i - 1] * commonRatio;
    }
    return sequence;
}

#endif
