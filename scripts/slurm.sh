#!/bin/bash

#SBATCH --account=fastdensekernel
#SBATCH --job-name=test_htool_benchmarks
#SBATCH --output=/mnt/beegfs/workdir/virgile.dubos/htool_benchmarks/output/test.out
#SBATCH --error=/mnt/beegfs/workdir/virgile.dubos/htool_benchmarks/output/test.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=cpu_shared

# mail alert at start, end and abortion of execution
## load modules
module purge
module load cmake/3.26.3 
module load git/2.31.1 
module load gcc/10.2.0 
module load lapack/3.10.0  
module load openmpi/4.1.4  
module load openblas/0.3.15

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## exe
./../build/bench_hmatrix_build_vs_pbl_size