#include "domain.h"

#ifdef USE_MPI

#include <stdlib.h>

/* Las claves Morton usan 21 bits por eje => 63 bits: el espacio es [0, 2^63). */
#define KEY_LIMIT (1ULL << 63)

void domain_init(Domain *d, MPI_Comm comm) {
    MPI_Comm_rank(comm, &d->rank);
    MPI_Comm_size(comm, &d->nprocs);
    d->splitters = malloc((d->nprocs + 1) * sizeof(uint64_t));
    d->splitters[0] = 0;
    d->splitters[d->nprocs] = UINT64_MAX;
}

void domain_free(Domain *d) {
    free(d->splitters);
    d->splitters = NULL;
}

void domain_global_bounds(Domain *d, const Particle *p, int n_local, MPI_Comm comm) {
    double lmin[3] = { 1e30,  1e30,  1e30};
    double lmax[3] = {-1e30, -1e30, -1e30};

    /* Un proceso sin partículas aporta los neutros de MIN/MAX y no altera el
       resultado, que es lo que hace seguro el caso n_local == 0. */
    for (int i = 0; i < n_local; i++) {
        for (int k = 0; k < 3; k++) {
            if (p[i].pos[k] < lmin[k]) lmin[k] = p[i].pos[k];
            if (p[i].pos[k] > lmax[k]) lmax[k] = p[i].pos[k];
        }
    }

    double gmin[3], gmax[3];
    MPI_Allreduce(lmin, gmin, 3, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(lmax, gmax, 3, MPI_DOUBLE, MPI_MAX, comm);

    /* Cubo con margen del 1%, misma construcción que octree_bounds para que el
       árbol y el particionamiento vean la misma geometría. */
    double side = 0.0;
    for (int k = 0; k < 3; k++)
        if (gmax[k] - gmin[k] > side) side = gmax[k] - gmin[k];
    if (side <= 0.0) side = 1.0;
    side += 2.0 * (side * 0.01);

    for (int k = 0; k < 3; k++) {
        double c = 0.5 * (gmin[k] + gmax[k]);
        d->gmin[k] = c - 0.5 * side;
        d->gmax[k] = c + 0.5 * side;
    }
}

/*
 * Cantidad de partículas locales con clave estrictamente menor que key.
 * Búsqueda binaria: asume el arreglo local ordenado por clave creciente.
 */
static int64_t count_below(const Particle *p, int n, uint64_t key) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (p[mid].morton < key) lo = mid + 1;
        else                     hi = mid;
    }
    return (int64_t)lo;
}

void domain_partition(Domain *d, const Particle *p, int n_local,
                      int64_t n_global, MPI_Comm comm) {
    domain_partition_ex(d, p, n_local, n_global, 0, comm);
}

void domain_partition_ex(Domain *d, const Particle *p, int n_local,
                         int64_t n_global, int use_work, MPI_Comm comm) {
    const int P = d->nprocs;
    d->splitters[0] = 0;
    d->splitters[P] = UINT64_MAX;
    if (P == 1) return;

    const int ns = P - 1;   /* splitters interiores a determinar */

    /*
     * Prefijo de pesos: prefix[j] = trabajo acumulado de las j primeras
     * partículas locales. Como el arreglo está ordenado por clave, el trabajo
     * con clave < K es prefix[count_below(K)] — misma búsqueda binaria que el
     * conteo, con una lectura extra. Sin pesos el prefijo es la identidad y se
     * usa directamente el índice.
     */
    int64_t *prefix = NULL;
    int64_t  total  = n_global;

    if (use_work) {
        prefix = malloc((size_t)(n_local + 1) * sizeof(int64_t));
        prefix[0] = 0;
        for (int i = 0; i < n_local; i++)
            prefix[i + 1] = prefix[i] + p[i].work;

        int64_t w_local = prefix[n_local], w_global;
        MPI_Allreduce(&w_local, &w_global, 1, MPI_INT64_T, MPI_SUM, comm);

        /* Primer paso: work todavía no fue medido. Caer al reparto por conteo
           en vez de dividir por cero o dejar todos los splitters en 0. */
        if (w_global <= 0) {
            free(prefix);
            prefix = NULL;
            use_work = 0;
        } else {
            total = w_global;
        }
    }

    uint64_t *lo     = malloc(ns * sizeof(uint64_t));
    uint64_t *hi     = malloc(ns * sizeof(uint64_t));
    uint64_t *mid    = malloc(ns * sizeof(uint64_t));
    int64_t  *target = malloc(ns * sizeof(int64_t));
    int64_t  *cl     = malloc(ns * sizeof(int64_t));
    int64_t  *cg     = malloc(ns * sizeof(int64_t));

    /*
     * splitter[k+1] = la menor clave S tal que (# partículas con clave < S)
     * alcanza el objetivo (k+1)·N/P. Bisección sobre S en [0, 2^63).
     */
    for (int k = 0; k < ns; k++) {
        lo[k] = 0;
        hi[k] = KEY_LIMIT;
        target[k] = ((int64_t)(k + 1) * total) / P;
    }

    /*
     * Los ns splitters se bisectan a la vez: una sola Allreduce de ns enteros
     * por iteración en lugar de ns Allreduce independientes. Como lo[] y hi[]
     * se derivan solo de valores globales, todos los procesos recorren
     * exactamente las mismas iteraciones y la condición de corte es idéntica en
     * todos: no hay riesgo de que un proceso salga del bucle colectivo antes.
     */
    for (int it = 0; it < 63; it++) {
        int active = 0;
        for (int k = 0; k < ns; k++) {
            if (lo[k] < hi[k]) {
                mid[k] = lo[k] + (hi[k] - lo[k]) / 2;
                active = 1;
            } else {
                mid[k] = lo[k];
            }
            int64_t j = count_below(p, n_local, mid[k]);
            cl[k] = prefix ? prefix[j] : j;
        }
        if (!active) break;

        MPI_Allreduce(cl, cg, ns, MPI_INT64_T, MPI_SUM, comm);

        for (int k = 0; k < ns; k++) {
            if (lo[k] >= hi[k]) continue;
            if (cg[k] >= target[k]) hi[k] = mid[k];
            else                    lo[k] = mid[k] + 1;
        }
    }

    for (int k = 0; k < ns; k++)
        d->splitters[k + 1] = lo[k];

    free(lo); free(hi); free(mid); free(target); free(cl); free(cg);
    free(prefix);
}

int domain_owner(const Domain *d, uint64_t key) {
    /* Mayor k con splitters[k] <= key. splitters[0] == 0 siempre cumple, así
       que el resultado está bien definido para cualquier clave. */
    int res = 0, l = 1, r = d->nprocs;
    while (l < r) {
        int m = l + (r - l) / 2;
        if (d->splitters[m] <= key) { res = m; l = m + 1; }
        else                        { r = m; }
    }
    return res;
}

#endif /* USE_MPI */
