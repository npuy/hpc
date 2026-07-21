CC = gcc
MPICC = mpicc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS = -lm

# OpenMP: el gcc de Apple es clang y necesita libomp de homebrew con -Xpreprocessor.
# En Linux (maquinas pcunix de fing) alcanza con -fopenmp.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(LIBOMP),)
    OMP_CFLAGS  = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    OMP_LDFLAGS = -L$(LIBOMP)/lib -lomp
  endif
else
  OMP_CFLAGS  = -fopenmp
  OMP_LDFLAGS = -fopenmp
endif

SRC = src/main.c src/particle.c src/force.c src/leapfrog.c src/morton.c \
      src/octree.c src/barnes_hut.c src/validation.c
HDR = src/vec3.h src/particle.h src/force.h src/leapfrog.h src/morton.h \
      src/octree.h src/barnes_hut.h src/validation.h

# Modulos que solo existen en la version MPI (semana 3).
MPI_SRC = src/mpi_types.c src/domain.c src/migration.c src/metrics.c src/let.c
MPI_HDR = src/mpi_types.h src/domain.h src/migration.h src/metrics.h src/let.h

.PHONY: all clean

all: nbody nbody_mpi

# Version secuencial + OpenMP (sin MPI).
nbody: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(OMP_CFLAGS) -o $@ $(SRC) $(LDFLAGS) $(OMP_LDFLAGS)

# Version hibrida MPI + OpenMP. USE_MPI activa las guardas del codigo distribuido.
nbody_mpi: $(SRC) $(HDR) $(MPI_SRC) $(MPI_HDR)
	$(MPICC) $(CFLAGS) $(OMP_CFLAGS) -DUSE_MPI -o $@ $(SRC) $(MPI_SRC) $(LDFLAGS) $(OMP_LDFLAGS)

clean:
	rm -f nbody nbody_mpi *.csv
