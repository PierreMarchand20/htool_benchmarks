#ifndef HTOOL_BENCHMARKS_MISC_HPP
#define HTOOL_BENCHMARKS_MISC_HPP

#include "xsimd/xsimd.hpp"
#include <complex>
#include <string>

template <typename T, int vectorisation>
struct allocator_type_spec {
    typedef typename std::allocator<T> allocator;
};

template <typename T>
struct allocator_type_spec<T, 3> {
    typedef xsimd::aligned_allocator<T> allocator;
};

template <class T, int vectorisation>
using allocator_type = typename allocator_type_spec<T, vectorisation>::allocator;

std::string vectorisation_name(int i) {
    std::string vec_name;
    switch (i) {
    case 0:
        vec_name = "none";
        break;
    case 1:
        vec_name = "vcl";
        break;
    case 2:
        vec_name = "xsimd (unaligned)";
        break;
    case 3:
        vec_name = "xsimd (aligned)";
        break;
    default:
        break;
    }
    return vec_name;
}

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

#endif