#!/bin/bash

# Initialization
MY_PATH="`dirname \"$0\"`"             
MY_PATH="`( cd \"$MY_PATH\" && pwd )`" 
cd ${MY_PATH}
cd ..

# Build
mkdir -p build
cd build
rm bench_hmatrix.csv
CC=gcc CXX=g++ cmake ../ -DCMAKE_BUILD_TYPE=Release_w_vecto
make BenchHMatrix
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
executable=${buildpath}/BenchHMatrix
time=00:30:00

# Log folder
mkdir -p ${MY_PATH}/../logs
logpath=${MY_PATH}/../logs

# Scripts
sizes=(1000000)
types=("S" "D")
vectorisations=(0 2 3)
compressors=("sympartialACA")
clusterings=("PCARegularClustering")

for ntask in "${ntasks[@]}"
do
  for thread in  "${threads[@]}"
  do
    for T in "${types[@]}"
    do
      for clustering in "${clusterings[@]}"
      do  
        for size in "${sizes[@]}"
        do
          for vectorisation in "${vectorisations[@]}"
          do
            for compressor in "${compressors[@]}"
            do  
                if ! ([ $T == "S" -a ${vectorisation} -eq 2 ] ||  [ $T == "S" -a ${vectorisation} -eq 3 ])
                then
                  if ((${ntask}*${thread}<=${nprocs}))
                  then 
                    nnodes=$((1+${ntask}*${thread}/${procs_per_node}))
                    if ((${ntask}*${thread} % ${procs_per_node}==0))
                    then
                      nnodes=$((${nnodes}-1))
                    fi
                    signature=${nnodes}_${ntask}_${thread}_${T}_${clustering}_${size}_${vectorisation}_${compressor}
                    echo ${MY_PATH}/slurm_bench_hmatrix.sh ${nnodes} ${ntask} $((ntask/nnodes)) ${thread} ${time} ${signature} ${executable} ${T} ${clustering} ${size} ${vectorisation} ${compressor} ${logpath} ${outputpath}
                    ${MY_PATH}/slurm_bench_hmatrix.sh ${nnodes} ${ntask} $((ntask/nnodes)) ${thread} ${time} ${signature} ${executable} ${T} ${clustering} ${size} ${vectorisation} ${compressor} ${logpath} ${outputpath}
                  fi
                fi
            done
          done
        done
      done
    done
  done
done

complex_types=("C" "Z")
sizes=(100000)
for ntask in "${ntasks[@]}"
do
  for thread in  "${threads[@]}"
  do
    for T in "${complex_types[@]}"
    do
      for clustering in "${clusterings[@]}"
      do  
        for size in "${sizes[@]}"
        do
          for vectorisation in "${vectorisations[@]}"
          do
            for compressor in "${compressors[@]}"
            do  
                if ! ([ $T == "C" -a $vectorisation -eq 2 ] ||  [ $T == "C" -a $vectorisation -eq 3 ])
                then
                  if ((${ntask}*${thread}<=${nprocs}))
                  then 
                    nnodes=$((1+${ntask}*${thread}/${procs_per_node}))
                    if ((${ntask}*${thread} % ${procs_per_node}==0))
                    then
                      nnodes=$((${nnodes}-1))
                    fi
                    signature=${nnodes}_${ntask}_${thread}
                    echo ${MY_PATH}/slurm_bench_hmatrix.sh ${nnodes} ${ntask} $((ntask/nnodes)) ${thread} ${time} ${signature} ${executable} ${T} ${clustering} ${size} ${vectorisation} ${compressor} ${logpath} ${outputpath}
                    ${MY_PATH}/slurm_bench_hmatrix.sh ${nnodes} ${ntask} $((ntask/nnodes)) ${thread} ${time} ${signature} ${executable} ${T} ${clustering} ${size} ${vectorisation} ${compressor} ${logpath} ${outputpath}
                  fi
                fi
            done
          done
        done
      done
    done
  done
done
