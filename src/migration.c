#include "migration.h"

#ifdef USE_MPI

#include "morton.h"
#include <stdlib.h>
#include <string.h>

int migrate_particles(Particle **p, int n_local, int *capacity,
                      const Domain *d, MPI_Datatype ptype, MPI_Comm comm) {
    const int P = d->nprocs;

    /* Copias no-const: morton_encode recibe los límites por puntero mutable. */
    double gmin[3], gmax[3];
    for (int k = 0; k < 3; k++) { gmin[k] = d->gmin[k]; gmax[k] = d->gmax[k]; }

    /* 1. Las partículas se movieron en el drift: recodificar sus claves con el
          bounding box global congelado. */
    morton_keys(*p, n_local, gmin, gmax);

    /* 2-3. Contar cuántas van a cada destino y calcular los desplazamientos. */
    int *dest    = malloc((n_local > 0 ? n_local : 1) * sizeof(int));
    int *scount  = calloc(P, sizeof(int));
    int *sdispl  = malloc(P * sizeof(int));
    int *rcount  = malloc(P * sizeof(int));
    int *rdispl  = malloc(P * sizeof(int));
    int *cursor  = malloc(P * sizeof(int));

    for (int i = 0; i < n_local; i++) {
        dest[i] = domain_owner(d, (*p)[i].morton);
        scount[dest[i]]++;
    }
    sdispl[0] = 0;
    for (int k = 1; k < P; k++) sdispl[k] = sdispl[k-1] + scount[k-1];

    /* 4. Agrupar por destino (counting sort estable sobre dest[]). */
    memcpy(cursor, sdispl, P * sizeof(int));
    Particle *sendbuf = malloc((n_local > 0 ? n_local : 1) * sizeof(Particle));
    for (int i = 0; i < n_local; i++)
        sendbuf[cursor[dest[i]]++] = (*p)[i];

    /* 5. Cada proceso averigua cuánto va a recibir de cada otro. */
    MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, comm);

    rdispl[0] = 0;
    for (int k = 1; k < P; k++) rdispl[k] = rdispl[k-1] + rcount[k-1];
    int n_recv = rdispl[P-1] + rcount[P-1];

    /* 6. Intercambio. n_recv puede ser 0: es legal y Alltoallv lo maneja. */
    Particle *recvbuf = malloc((n_recv > 0 ? n_recv : 1) * sizeof(Particle));
    MPI_Alltoallv(sendbuf, scount, sdispl, ptype,
                  recvbuf, rcount, rdispl, ptype, comm);

    free(*p);
    *p = recvbuf;
    *capacity = (n_recv > 0 ? n_recv : 1);

    /* 7. Restaurar la invariante de orden por clave. Recodificar con el mismo
          box congelado es idempotente (las posiciones no cambiaron en el
          intercambio), así que esto solo reordena. */
    morton_sort(*p, n_recv, gmin, gmax);

    free(sendbuf); free(dest); free(scount);
    free(sdispl); free(rcount); free(rdispl); free(cursor);

    return n_recv;
}

void replicate_particles(const Particle *local, int n_local,
                         Particle *all, int n_global, int *offset,
                         MPI_Datatype ptype, MPI_Comm comm) {
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));

    MPI_Allgather(&n_local, 1, MPI_INT, counts, 1, MPI_INT, comm);

    displs[0] = 0;
    for (int k = 1; k < P; k++) displs[k] = displs[k-1] + counts[k-1];
    *offset = displs[rank];

    (void)n_global;  /* el caller garantiza capacidad; Σ counts == n_global */
    MPI_Allgatherv(local, n_local, ptype, all, counts, displs, ptype, comm);

    free(counts);
    free(displs);
}

void particle_checksum(const Particle *p, int n_local,
                       int64_t *sum_id, int64_t *sum_id2, MPI_Comm comm) {
    int64_t local[2] = {0, 0};
    for (int i = 0; i < n_local; i++) {
        local[0] += p[i].id;
        local[1] += p[i].id * p[i].id;
    }
    int64_t global[2];
    MPI_Allreduce(local, global, 2, MPI_INT64_T, MPI_SUM, comm);
    *sum_id  = global[0];
    *sum_id2 = global[1];
}

#endif /* USE_MPI */
