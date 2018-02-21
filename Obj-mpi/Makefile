# -*- Makefile -*-
SHELL=/bin/sh
CC=mpicxx
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -lnetcdf -std=c++11 -fopenmp
LDLIBS=-lm -lnetcdf

# list of source files
SRC=ljmd-mpi.c timers.c ncdf.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-mpi.x

clean:
	rm -f *.mod *.o ../ljmd-mpi.x

# linker rule
../ljmd-mpi.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
ncdf.o: ncdf.c ncdf.h
timers.o: timers.c timers.h
ljmd.o: ljmd.c timers.h
