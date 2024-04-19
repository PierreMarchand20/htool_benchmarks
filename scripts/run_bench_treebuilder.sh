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

# change CPU mode to performance BEFORE running benchmarks
sudo cpupower frequency-set --governor performance > /dev/null

# Scripts
./bench_hmatrix_build

# change back CPU mode to powerersave AFTER running benchmarks
sudo cpupower frequency-set --governor powersave > /dev/null

# check current CPU mode
# cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor