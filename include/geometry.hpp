#ifndef HTOOL_BENCHMARKS_GEOMETRY_HPP
#define HTOOL_BENCHMARKS_GEOMETRY_HPP

#include <array>
#include <cmath>

double create_geometry(int size, double *points, std::array<double, 3> shift = {0, 0, 0}) {
    double radius, step, k;
    radius                 = 1.;
    step                   = 1.75 * M_PI * radius / sqrt((double)size);
    k                      = 2 * M_PI / (10. * step); // 10 points / lambda
    double length          = 2 * M_PI * radius;
    double pointsPerCircle = length / step;
    double angleStep       = 2 * M_PI / pointsPerCircle;
    for (int j = 0; j < size; j++) {
        points[3 * j + 0] = shift[0] + radius * cos(angleStep * j);
        points[3 * j + 1] = shift[1] + radius * sin(angleStep * j);
        points[3 * j + 2] = shift[2] + (step * j) / pointsPerCircle;
    }
    return k;
}

#endif