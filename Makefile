# LIBRARY_DIRS = /usr/lib/x86_64-linux-gnu/hdf5/serial
# INCLUDE_DIRS = /usr/include/hdf5/serial
LIBRARY_DIRS = ~/usr/lib
INCLUDE_DIRS = ~/usr/include
CFLAGS = -I$(INCLUDE_DIRS) -L$(LIBRARY_DIRS) -std=c++0x -Wall -lm -lhdf5 -O0

output:
	g++ sheep_v4.cpp ${CFLAGS} -o sheep
#sheep.o: sheep_v4.cpp animal.h
#	g++ -c ${CFLAGS} sheep_v4.cpp

#animal.o: animal.cpp animal.h
#    g++ -c ${CFLAGS} animal.cpp
clean:
	rm sheep 
