#ifndef HTOOL_BENCHMARKS_GEOMETRY
#define HTOOL_BENCHMARKS_GEOMETRY

#include <array>
#include <fstream>
#include <iostream>
#include <vector>

std::vector<std::array<double, 3>> LoadGeometry(std::string filename) {

    std::ifstream file(filename);

    if (!file.is_open())
        throw std::runtime_error("Could not open file");
    else {
        std::vector<std::array<double, 3>> geometry;
        std::string line;
        double coord;
        std::array<double, 3> point;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            for (int i = 0; i < 3; i++) {
                iss >> coord;
                point[i] = coord;
            }
            geometry.emplace_back(point);
        }

        return geometry;
    }
}
#endif