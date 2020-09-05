# LIBRARY_DIRS = /usr/lib/x86_64-linux-gnu/hdf5/serial
# INCLUDE_DIRS = /usr/include/hdf5/serial
LIBRARY_DIRS = ~/usr/lib
INCLUDE_DIRS = ~/usr/include

output:
	gcc sheep_v3.c -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -lhdf5 -std=gnu99 -lm -o sheep -O0
clean:
	rm sheep