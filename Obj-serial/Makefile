# -*- Makefile -*-
SHELL=/bin/sh
CC=g++
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=ljmd.c timers.h
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x

clean:
	rm -f *.mod *.o ../ljmd-serial.x

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<
%.gch: %.h
	$(CC) -c $(CFLAGS) $<


# dependencies
timers.gch: timers.h
ljmd.o: ljmd.c timers.h
