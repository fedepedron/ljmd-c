# -*- Makefile -*-
SHELL=/bin/sh
CC=mpicxx
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -lnetcdf -std=c++11 -fopenmp
LDLIBS=-lm -lnetcdf

# list of source files
SRC=ljmd-mpi.c timers.c md_utils.c md_io.c
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

timers.o: timers.h timers.c
md_utils.o: md_utils.h md_utils.c
md_io.o: md_io.h md_utils.h md_io.c
ljmd.o: timers.h md_io.h md_utils.h ljmd.c
