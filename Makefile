LIBRARY_DIRS = ~/usr/lib
INCLUDE_DIRS = ~/usr/include
# CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS)  -lm -O0
CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -std=c++0x -Wall -lm -O0
# CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -lhdf5 -std=c++0x -Wall -lm -O0

output:
	g++-10 sheep_v4_omp.cpp ${CFLAGS} -o sheep
omp:
	g++-10 -fopenmp sheep_v4_omp.cpp ${CFLAGS} -o sheep
clean:
	rm sheep 
