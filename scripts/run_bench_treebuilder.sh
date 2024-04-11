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
make BenchTreeBuilder
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
executable=${buildpath}/BenchTreeBuilder
time=00:30:00

# Log folder
mkdir -p ${MY_PATH}/../logs
logpath=${MY_PATH}/../logs

# Scripts
./BenchTreeBuilder