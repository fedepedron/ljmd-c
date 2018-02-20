# -*- Makefile -*-
SHELL=/bin/sh
CC=mpic++
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fopenmp
LDLIBS=-lm

# list of source files
SRC=ljmd-mpi.c timers.h
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
%.gch: %.h
	$(CC) -c $(CFLAGS) $<


# dependencies
timers.gch: timers.h
ljmd.o: ljmd.c timers.h
