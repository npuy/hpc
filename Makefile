CC = gcc
MPICC = mpicc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS = -lm

SRC = src/main.c src/particle.c src/force.c src/leapfrog.c src/morton.c src/validation.c
HDR = src/vec3.h src/particle.h src/force.h src/leapfrog.h src/morton.h src/validation.h

.PHONY: all clean

all: nbody

nbody: $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

nbody_mpi: $(SRC) $(HDR)
	$(MPICC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

clean:
	rm -f nbody nbody_mpi *.csv
