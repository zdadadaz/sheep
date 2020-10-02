#!/bin/bash
#SBATCH −−partition=coursework
#SBATCH −−job−name=:sad_parrot:
#SBATCH −−nodes=1
#SBATCH −−ntasks=4
#SBATCH −−ntasks−per−node=4
#SBATCH −−cpus−per−task=4

export OMP_NUM_THREADS=4
export SLURM_TASKS_PER_NODE=4
export SLURM_NPROCS=4

DATE=$(date +"%Y%m%d%H%M")
echo "time started "$DATE
echo "This is job ’$SLURM_JOB_NAME’ (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "running with OMP_NUM_THREADS= $OMP_NUM_THREADS "
echo "running with SLURM_TASKS_PER_NODE= $SLURM_TASKS_PER_NODE "
echo "running with SLURM_NPROCS= $SLURM_NPROCS "
echo "Now we start the show:"
export TIMEFORMAT="%E sec"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NPROCS} ./bin/lcs-hybrid test_input.txt test_output.txt

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
# echo "we just ran with the following SLURM environment variables" # env | grep SLURM
