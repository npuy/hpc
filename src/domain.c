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

    d->orb       = malloc((size_t)(2 * d->nprocs) * sizeof(OrbNode));
    d->orb_count = 0;
    d->use_orb   = 0;
}

void domain_free(Domain *d) {
    free(d->splitters);
    free(d->orb);
    d->splitters = NULL;
    d->orb = NULL;
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

/* ------------------------------------------------------------------ */
/* ORB: biseccion recursiva ortogonal (semana 5)                       */
/* ------------------------------------------------------------------ */

/* Peso de una particula: su costo medido, o 1 si todavia no se midio. */
static int64_t weight_of(const Particle *p, int use_work) {
    return use_work ? p->work : 1;
}

/*
 * Construye recursivamente el subarbol que cubre los ranks
 * [first, first+nprocs) sobre la caja [bmin, bmax].
 *
 * idx[0..n_idx) son los indices de las particulas LOCALES que caen en esta caja.
 * La pertenencia se lleva por el camino recorrido y no por un test geometrico:
 * asi el conteo de la biseccion y domain_owner_pos usan exactamente el mismo
 * predicado, que es la condicion para no perder particulas sobre el plano.
 *
 * Todos los procesos recorren la misma estructura de recursion, porque depende
 * solo de valores globales. Es lo que permite usar colectivos aca adentro.
 */
static int orb_build(Domain *d, const Particle *p, int *idx, int n_idx,
                     const double bmin[3], const double bmax[3],
                     int first, int nprocs, int use_work, MPI_Comm comm) {
    int self = d->orb_count++;
    for (int k = 0; k < 3; k++) {
        d->orb[self].min[k] = bmin[k];
        d->orb[self].max[k] = bmax[k];
    }
    d->orb[self].first  = first;
    d->orb[self].nprocs = nprocs;

    if (nprocs == 1) {
        d->orb[self].axis  = -1;
        d->orb[self].left  = -1;
        d->orb[self].right = -1;
        return self;
    }

    /* Eje = dimension mas larga. Evita cajas con forma de lamina, que tienen
       mucha superficie y por lo tanto generan mucho LET. */
    int axis = 0;
    double ext = bmax[0] - bmin[0];
    for (int k = 1; k < 3; k++)
        if (bmax[k] - bmin[k] > ext) { ext = bmax[k] - bmin[k]; axis = k; }

    int64_t w_local = 0;
    for (int j = 0; j < n_idx; j++) w_local += weight_of(&p[idx[j]], use_work);
    int64_t w_total;
    MPI_Allreduce(&w_local, &w_total, 1, MPI_INT64_T, MPI_SUM, comm);

    /* Reparto proporcional: soporta P que no es potencia de 2 sin parches. */
    int nL = nprocs / 2;
    int64_t target = (w_total * nL) / nprocs;

    /* Biseccion sobre la coordenada. Corte por numero fijo de iteraciones y no
       por tolerancia: asi todos los procesos hacen exactamente las mismas
       llamadas colectivas. 50 iteraciones agotan la precision de un double. */
    double lo = bmin[axis], hi = bmax[axis], mid = 0.5 * (lo + hi);
    for (int it = 0; it < 50; it++) {
        mid = 0.5 * (lo + hi);
        int64_t below_local = 0;
        for (int j = 0; j < n_idx; j++)
            if (p[idx[j]].pos[axis] < mid)
                below_local += weight_of(&p[idx[j]], use_work);
        int64_t below;
        MPI_Allreduce(&below_local, &below, 1, MPI_INT64_T, MPI_SUM, comm);

        if (below < target) lo = mid;
        else                hi = mid;
    }
    double split = 0.5 * (lo + hi);

    d->orb[self].axis  = axis;
    d->orb[self].split = split;

    /* Particion in-place de idx[] por el mismo predicado (<) que usara
       domain_owner_pos. Cada nivel cuesta O(n_idx), no O(n_local) por nodo. */
    int nleft = 0;
    for (int j = 0; j < n_idx; j++) {
        if (p[idx[j]].pos[axis] < split) {
            int tmp = idx[nleft]; idx[nleft] = idx[j]; idx[j] = tmp;
            nleft++;
        }
    }

    double lmax[3], rmin[3];
    for (int k = 0; k < 3; k++) { lmax[k] = bmax[k]; rmin[k] = bmin[k]; }
    lmax[axis] = split;
    rmin[axis] = split;

    int l = orb_build(d, p, idx,          nleft,         bmin, lmax,
                      first,      nL,           use_work, comm);
    int r = orb_build(d, p, idx + nleft,  n_idx - nleft, rmin, bmax,
                      first + nL, nprocs - nL,  use_work, comm);

    d->orb[self].left  = l;
    d->orb[self].right = r;
    return self;
}

void domain_partition_orb(Domain *d, const Particle *p, int n_local,
                          int use_work, MPI_Comm comm) {
    d->orb_count = 0;
    d->use_orb   = 1;

    /* Primer paso: work vale 0 para todas. Caer al reparto por conteo en vez de
       bisecar un total nulo y dejar todos los cortes pegados a un borde. */
    if (use_work) {
        int64_t w = 0, w_global;
        for (int i = 0; i < n_local; i++) w += p[i].work;
        MPI_Allreduce(&w, &w_global, 1, MPI_INT64_T, MPI_SUM, comm);
        if (w_global <= 0) use_work = 0;
    }

    int *idx = malloc((size_t)(n_local > 0 ? n_local : 1) * sizeof(int));
    for (int i = 0; i < n_local; i++) idx[i] = i;

    orb_build(d, p, idx, n_local, d->gmin, d->gmax, 0, d->nprocs, use_work, comm);

    free(idx);
}

int domain_owner_pos(const Domain *d, const double pos[3]) {
    int i = 0;
    while (d->orb[i].axis >= 0)
        i = (pos[d->orb[i].axis] < d->orb[i].split) ? d->orb[i].left
                                                    : d->orb[i].right;
    return d->orb[i].first;
}

void domain_box(const Domain *d, int rank, double min[3], double max[3]) {
    for (int i = 0; i < d->orb_count; i++) {
        if (d->orb[i].axis < 0 && d->orb[i].first == rank) {
            for (int k = 0; k < 3; k++) {
                min[k] = d->orb[i].min[k];
                max[k] = d->orb[i].max[k];
            }
            return;
        }
    }
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
