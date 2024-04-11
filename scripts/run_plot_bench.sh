#!/bin/bash

# Initialization
MY_PATH="`dirname \"$0\"`"             
MY_PATH="`( cd \"$MY_PATH\" && pwd )`" 
cd ${MY_PATH}
cd ..

# Build
mkdir -p build
cd build

# Scripts
python ../src/plot_bench.py bench_hmatrix_build_vs_pbl_size.csv 1
python ../src/plot_bench.py bench_hmatrix_build_vs_thread.csv 1
python ../src/plot_bench.py bench_hmatrix_build_vs_ratio.csv 1
python ../src/plot_bench.py bench_hmatrix_matrix_product_vs_pbl_size.csv 1
python ../src/plot_bench.py bench_hmatrix_matrix_product_vs_thread.csv 1
python ../src/plot_bench.py bench_hmatrix_matrix_product_vs_ratio.csv 1
python ../src/plot_bench.py bench_hmatrix_factorization_LU_vs_pbl_size.csv 1
python ../src/plot_bench.py bench_hmatrix_factorization_Cholesky_vs_pbl_size.csv 1

