#!/bin/bash
# Initialization
MY_PATH="`dirname \"$0\"`"             
MY_PATH="`( cd \"$MY_PATH\" && pwd )`" 
cd ${MY_PATH}
cd ..

# Build
mkdir -p build
cd build
CC=gcc CXX=g++ cmake ../ -DCMAKE_BUILD_TYPE=Release_w_vecto
make

# Setup
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close 
export OMP_PLACES=cores


sizes=(1000000)
threads=(1 2 4 8 16)
compressors=("SVD" "fullACA" "partialACA" "sympartialACA")
clusterings=("PCARegularClustering" "PCAGeometricClustering" "BoundingBox1RegularClustering" "BoundingBox1GeometricClustering")

# BenchCompression
types=("S" "D" "C" "Z")
vectorisations=(0 1)
export OMP_NUM_THREADS=1
for T in "${types[@]}"
do  
    for clustering in "${clusterings[@]}"
    do  
        for vectorisation in "${vectorisations[@]}"
        do  
            for compressor in "${compressors[@]}"
            do  
                ./BenchCompression 1000 1000 ${clustering} ${T} ${vectorisation} ${compressor}
                ./BenchCompression 10000 2500 ${clustering} ${T} ${vectorisation} ${compressor}
            done
        done
    done
done

types=("D" "Z")
vectorisations=(2 3)
for T in "${types[@]}"
do  
    for clustering in "${clusterings[@]}"
    do  
        for vectorisation in "${vectorisations[@]}"
        do  
            for compressor in "${compressors[@]}"
            do  
                ./BenchCompression 1000 1000 ${clustering} ${T} ${vectorisation} ${compressor}
                ./BenchCompression 10000 2500 ${clustering} ${T} ${vectorisation} ${compressor}
            done
        done
    done
done


# Output
cd ..
mkdir -p output
cp build/*.csv output/

