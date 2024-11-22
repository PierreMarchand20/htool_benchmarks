#!/bin/bash

#SBATCH --account=fastdensekernel
#SBATCH --job-name=bench_hmatrix_build
#SBATCH --output=/mnt/beegfs/workdir/virgile.dubos/htool_benchmarks/output/bench_hmatrix_build.out
#SBATCH --error=/mnt/beegfs/workdir/virgile.dubos/htool_benchmarks/output/bench_hmatrix_build.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=cpu_shared

# mail alert at start, end and abortion of execution
## load modules
## module load cmake/3.26.3 git/2.31.1 gcc/10.2.0 lapack/3.10.0  openmpi/4.1.4  openblas/0.3.15
module purge
module load cmake/3.26.3 
module load git/2.31.1 
module load gcc/10.2.0 
module load lapack/3.10.0  
module load openmpi/4.1.4  
module load openblas/0.3.15

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## exe
./../build/bench_hmatrix_build 0
./../build/bench_hmatrix_build 1