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
# nprocs=16
# ntasks=(1 2 4 8 16)
# threads=(1 2 4 8 16)
# procs_per_node=16

# Misc datas
executable=${buildpath}/bench_hmatrix_build
time=00:30:00

# Log folder
mkdir -p ${MY_PATH}/../logs
logpath=${MY_PATH}/../logs

# Prepare machine for benchmarking 
sudo cpupower frequency-set --governor performance > /dev/null # Change all CPUs mode to performance
sudo bash -c "echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo" # Disable turbo
sudo bash -c "echo 0 > /proc/sys/kernel/randomize_va_space" # Disable ASLR

# Scripts
## Benchmark options : 
### "--benchmark_out=<filename> --benchmark_out_format={json|console|csv}" write benchmark results to a file in the setting format
### "--benchmark_filter=<regex>" only run the benchmarks that match regex, e.g. --benchmark_filter=bench_hmatrix_build/128 
### "--benchmark_min_warmup_time=<seconds>"  warming up all benchmarks for at least this many seconds, can be overwrited by MinWarmUpTime for specific benchmark
### "--benchmark_repetitions=<number>" number of repetitions for all benchmarks, can be overwrited by Repetitions for specific benchmark
### "--benchmark_enable_random_interleaving" lower run-to-run variance by randomly interleaving repetitions of a microbenchmark with repetitions from other microbenchmarks in the same benchmark test.
### "taskset -c 0 ./bench_hmatrix_build" Set the benchmark program's task affinity to a fixed cpu e.g 0.

taskset -c 0 ./bench_hmatrix_build --benchmark_out=bench_hmatrix_build.json --benchmark_out_format=json --benchmark_enable_random_interleaving=true
./../external/benchmark/tools/compare.py filters bench_hmatrix_build.json BM_test_hmatrix_build BM_NEW_test_hmatrix_build # Compare two different filters of one benchmark

# Restore machine settings
sudo cpupower frequency-set --governor powersave > /dev/null # Change back CPU mode to powersave 
sudo bash -c "echo 0 > /sys/devices/system/cpu/intel_pstate/no_turbo" # Enable turbo
sudo bash -c "echo 2 > /proc/sys/kernel/randomize_va_space" # Enable ASLR

# Check current CPU setting
#  cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ; cat /sys/devices/system/cpu/intel_pstate/no_turbo # 0/1 == turbo enabled/disabled
#  cat /proc/sys/kernel/randomize_va_space # 0/1/2 == No randomization/Conservative randomization/Full randomization

# for i in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
#     do
#         cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
#     done
    
