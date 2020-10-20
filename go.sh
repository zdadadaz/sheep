#!/bin/bash
#SBATCH --job-name=jc-
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G  
#SBATCH --time=0-1:00 
#SBATCH --partition=cosc
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export TIMEFORMAT="%E sec"
echo "n ${SLURM_NNODES} tpn ${SLURM_TASKS_PER_NODE} t ${SLURM_NTASKS} ct ${SLURM_CPUS_PER_TASK}"
module load gnu
module load mpi/openmpi-x86_64

time mpirun -mca btl ^openib ./sheep_mpi

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
