#!/bin/bash
#SBATCH --job-name=jc-
#SBATCH --nodes=
#SBATCH --ntasks-per-node=
#SBATCH --ntasks=
#SBATCH --cpus-per-task=
#SBATCH --mem-per-cpu=2G  
#SBATCH --time=0-1:00 
#SBATCH --partition=cosc
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export TIMEFORMAT="%E sec"
echo "n ${SLURM_NNODES} tpn ${SLURM_TASKS_PER_NODE} t ${SLURM_NTASKS} ct ${SLURM_CPUS_PER_TASK}"
module load gnu
module load mpi/openmpi-x86_64

time mpirun -mca btl ^openib ./tut07

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
