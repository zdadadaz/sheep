LIBRARY_DIRS = ~/usr/lib
INCLUDE_DIRS = ~/usr/include
# CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS)  -lm -O0
CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -std=c++0x -Wall -lm -O0
# CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -lhdf5 -std=c++0x -Wall -lm -O0

all: serial omp mpi hyb

serial:
	g++ sheep_v5.1_serial.cpp ${CFLAGS} -o sheep
omp:
	g++ -fopenmp sheep_v6.1_omp.cpp ${CFLAGS} -o sheep_omp
mpi:
	mpic++ sheep_v6.1_mpi.cpp ${CFLAGS} -o sheep_mpi
hyb:
	mpic++ -fopenmp  sheep_v6.1_hyb.cpp ${CFLAGS} -o sheep_hyb
clean:
	rm sheep sheep_omp sheep_mpi sheep_hyb

# time mpirun -n 1 --mca btl ^openib sheep_mpi
# module load mpi/openmpi-x86_64
# module load gnu
# squeue -u s4575321
# export OMP_NUM_THREADS=4