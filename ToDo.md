This file lists future features waiting to be implemented, in no particular order of priority.

- refactor the 3 python display scripts: plot_bench_vs_*.py to make them more robust. In particular, allow the user to plot graphs corresponding to a subset of the parameter sets in the csv file.
- avoid duplicating the code in bench_hmatrix_build.hpp, allowing TreeBuilder Classic and TaskBased to be managed in the same way as in hmatrix_fixture.hpp
- clean up the directory by removing obsolete files (GoogleBenchmark and old script structure) as well as the GoogleBenchmark sub-module
- pass symmetry_type as an argument at runtime for the 3 benchmark types

