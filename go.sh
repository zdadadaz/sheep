#!/bin/bash
#SBATCH --job-name=jc
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G  # memory
#SBATCH --time=0-1:00       # time (D-HH:MM)
#SBATCH --partition=cosc

# export OMP_NUM_THREADS=2

export TIMEFORMAT="%E sec"
echo "n ${SLURM_NNODES} tpn ${SLURM_TASKS_PER_NODE} t ${SLURM_NTASKS} ct ${SLURM_CPUS_PER_TASK}"
module load gnu
module load mpi/openmpi-x86_64
# make mpi
time mpirun -mca btl ^openib ./tut07
time mpirun -mca btl ^openib ./sheep_mpi
# time mpirun -n ${SLURM_NTASKS} -npernode ${SLURM_TASKS_PER_NODE}  -mca btl ^openib ./tut07

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE