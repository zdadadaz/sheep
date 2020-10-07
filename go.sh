#!/bin/bash
#SBATCH −−partition=cosc
#SBATCH −−job−name=jc
#SBATCH −−nodes=1
#SBATCH −−ntasks−per−node=12
#SBATCH −−ntasks=12
#SBATCH −−cpus−per−task=8
#SBATCH --time=1:00:00
#SBATCH --mem=120000

export SLURM_NNODES=1
export SLURM_NTASKS=12
export SLURM_TASKS_PER_NODE=12
export SLURM_CPUS_PER_TASK=8
export OMP_NUM_THREADS=8
export TIMEFORMAT="%E sec"
echo "n ${SLURM_NNODES} tpn ${SLURM_NTASKS} t ${SLURM_TASKS_PER_NODE} ct ${SLURM_CPUS_PER_TASK}"
# module load gnu
module load mpi/openmpi-x86_64
make mpi

time mpirun -n ${SLURM_NTASKS} -mca btl ^openib ./sheep_mpi

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
