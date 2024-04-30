#!/bin/bash

# Initialization
MY_PATH="`dirname \"$0\"`"             
MY_PATH="`( cd \"$MY_PATH\" && pwd )`" 
cd ${MY_PATH}
cd ..

# Build
mkdir -p build
cd build
# rm bench_hmatrix.csv
CC=gcc CXX=g++ cmake ../ -DCMAKE_BUILD_TYPE=Release_w_vecto
make bench_hmatrix_build
buildpath=${MY_PATH}/../build_intel

# Output
mkdir -p ${MY_PATH}/../output
outputpath=${MY_PATH}/../output/

# HPC data
nprocs=16
ntasks=(1 2 4 8 16)
threads=(1 2 4 8 16)
procs_per_node=16

# Misc datas
executable=${buildpath}/bench_hmatrix_build
time=00:30:00

# Log folder
mkdir -p ${MY_PATH}/../logs
logpath=${MY_PATH}/../logs

# Change CPU mode to performance BEFORE running benchmarks
sudo cpupower frequency-set --governor performance > /dev/null

# Scripts
## Compiler options : 
### "--benchmark_out=<filename> --benchmark_out_format={json|console|csv}" write benchmark results to a file in the setting format
### "--benchmark_filter=<regex>" only run the benchmarks that match regex, e.g. --benchmark_filter=bench_hmatrix_build/128 
./bench_hmatrix_build --benchmark_out=bench_hmatrix_build.json --benchmark_out_format=json

# Change back CPU mode to powersave AFTER running benchmarks
sudo cpupower frequency-set --governor powersave > /dev/null

# Check current CPU mode
# cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor