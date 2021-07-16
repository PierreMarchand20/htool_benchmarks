
# Initialization
MY_PATH="`dirname \"$0\"`"             
MY_PATH="`( cd \"$MY_PATH\" && pwd )`" 
cd ${MY_PATH}
cd ..

# Build
mkdir -p build
cd build
CC=gcc-9 CXX=g++-9 cmake ../ -DCMAKE_BUILD_TYPE=Release_w_vecto
make

# Setup
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close 
export OMP_PLACES=cores

types=(S D Z D)
sizes=(1000000)
threads=(1 2 4 8 16)
procs=(1 2 4 8 16)
compressors=("SVD" "fullACA" "partialACA" "sympartialACA")
clusterings=("PCARegularClustering" "PCAGeometricClustering" "BoundingBox1RegularClustering" "BoundingBox1GeometricClustering")
vectorisations=(0 1 2 3)

# BenchCompression
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

# BenchHmatrix
sizes=(1000000)
for T in "${types[@]}"
do
    for thread in "${threads[@]}"
    do
        for size in "${sizes[@]}"
        do
            for proc in "${procs[@]}"
            do
                for vectorisation in "${vectorisations[@]}"
                do
                    for compressor in "${compressors[@]}"
                    do  
                        export OMP_NUM_THREADS=${thread}
                        mpirun -np ${proc} BenchHmatrix ${size} ${clustering} ${T} ${vectorisation} ${compressor}
                    done
                done
            done
        done
    done
done

complex_types=(C Z)
sizes=(100000)
for T in "${complex_types[@]}"
do
    for thread in "${threads[@]}"
    do
        for size in "${sizes[@]}"
        do
            for proc in "${procs[@]}"
            do
                for vectorisation in "${vectorisations[@]}"
                do
                    for compressor in "${compressors[@]}"
                    do  
                        export OMP_NUM_THREADS=${thread}
                        mpirun -np ${proc} BenchHmatrix ${size} ${clustering} ${T} ${vectorisation} ${compressor}
                    done
                done
            done
        done
    done
done

# Output
cd ..
mkdir -p output
cp build/*.csv output/

