#include "let.h"

#ifdef USE_MPI

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vec3.h"

/* Lista de partículas que crece por duplicación (buffer de exportación). */
typedef struct {
    Particle *buf;
    int       n;
    int       cap;
} PList;

static void plist_init(PList *l) {
    l->cap = 256;
    l->n   = 0;
    l->buf = malloc((size_t)l->cap * sizeof(Particle));
}

static void plist_push(PList *l, const Particle *p) {
    if (l->n == l->cap) {
        l->cap *= 2;
        l->buf = realloc(l->buf, (size_t)l->cap * sizeof(Particle));
    }
    l->buf[l->n++] = *p;
}

void let_gather_domain_boxes(const Particle *local, int n_local,
                             const Domain *d, double *boxes, MPI_Comm comm) {
    double mine[6];

    if (n_local > 0) {
        for (int k = 0; k < 3; k++) {
            mine[k]     =  1e30;
            mine[3 + k] = -1e30;
        }
        for (int i = 0; i < n_local; i++) {
            for (int k = 0; k < 3; k++) {
                if (local[i].pos[k] < mine[k])     mine[k]     = local[i].pos[k];
                if (local[i].pos[k] > mine[3 + k]) mine[3 + k] = local[i].pos[k];
            }
        }
    } else {
        /* Caja degenerada en el centro del dominio: este proceso no calcula
           ninguna fuerza, así que no necesita fantasmas. Lo que le manden los
           demás es trabajo desperdiciado pero inofensivo. */
        for (int k = 0; k < 3; k++) {
            double c = 0.5 * (d->gmin[k] + d->gmax[k]);
            mine[k] = mine[3 + k] = c;
        }
    }

    MPI_Allgather(mine, 6, MPI_DOUBLE, boxes, 6, MPI_DOUBLE, comm);
}

/*
 * Distancia² mínima del punto v a la caja [bmin, bmax]. Cero si v cae dentro.
 * Es la cota que hace conservador al criterio de exportación: cualquier
 * partícula del destino está al menos a esta distancia del centro de masa.
 */
static double dist2_to_box(const double v[3],
                           const double bmin[3], const double bmax[3]) {
    double d2 = 0.0;
    for (int k = 0; k < 3; k++) {
        double dd = 0.0;
        if      (v[k] < bmin[k]) dd = bmin[k] - v[k];
        else if (v[k] > bmax[k]) dd = v[k] - bmax[k];
        d2 += dd * dd;
    }
    return d2;
}

/*
 * Recorre el subárbol idx decidiendo qué mandar al destino con caja
 * [bmin, bmax]: resumen (pseudo-partícula) si el nodo está lo bastante lejos,
 * descenso si no, y partículas reales al llegar a una hoja.
 */
static void select_nodes(const Octree *t, int idx, const Particle *p,
                         const double bmin[3], const double bmax[3],
                         double theta2, double eps2, PList *out) {
    const OctreeNode *nd = &t->nodes[idx];
    if (nd->mass == 0.0) return;   /* celda vacía: nada que exportar */

    if (nd->is_leaf) {
        for (int k = nd->particle; k != -1; k = t->next[k])
            plist_push(out, &p[k]);
        return;
    }

    double d2 = dist2_to_box(nd->cm, bmin, bmax);
    double s  = nd->max[0] - nd->min[0];

    /* Mismo MAC que accumulate_force, pero con la distancia mínima a la caja en
       lugar de la distancia a una partícula concreta. Con theta = 0 nunca se
       cumple, se desciende siempre hasta las hojas y el LET degenera en réplica
       completa: es lo que hace exacto al test 14. */
    if (s * s < theta2 * (d2 + eps2)) {
        Particle g;
        memset(&g, 0, sizeof(Particle));
        VEC3_COPY(g.pos, nd->cm);
        g.mass = nd->mass;
        g.id   = PARTICLE_GHOST_ID;
        plist_push(out, &g);
        return;
    }

    for (int ch = 0; ch < 8; ch++) {
        int c = nd->children[ch];
        if (c != -1)
            select_nodes(t, c, p, bmin, bmax, theta2, eps2, out);
    }
}

int let_exchange(Particle **p, int n_local, int *capacity,
                 const Octree *local_tree, const Domain *d,
                 double theta, double softening,
                 MPI_Datatype ptype, MPI_Comm comm) {
    const int P    = d->nprocs;
    const int rank = d->rank;

    double *boxes = malloc((size_t)6 * P * sizeof(double));
    let_gather_domain_boxes(*p, n_local, d, boxes, comm);

    double theta2 = theta * theta;
    double eps2   = softening * softening;

    int *scount = calloc(P, sizeof(int));
    int *sdispl = malloc(P * sizeof(int));
    int *rcount = malloc(P * sizeof(int));
    int *rdispl = malloc(P * sizeof(int));

    /* Los destinos se recorren en orden, así que el buffer queda ya agrupado
       por destino y los desplazamientos de envío salen de los conteos. */
    PList send;
    plist_init(&send);

    for (int r = 0; r < P; r++) {
        sdispl[r] = send.n;
        if (r == rank) continue;    /* nadie se importa a sí mismo */
        select_nodes(local_tree, local_tree->root, *p,
                     &boxes[6 * r], &boxes[6 * r + 3], theta2, eps2, &send);
        scount[r] = send.n - sdispl[r];
    }

    MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, comm);

    rdispl[0] = 0;
    for (int k = 1; k < P; k++) rdispl[k] = rdispl[k-1] + rcount[k-1];
    int n_ghost = rdispl[P-1] + rcount[P-1];

    /* Los fantasmas van después de las locales, en el mismo arreglo: así el
       árbol fusionado se arma sin copiar nada y compute_forces_bh_range calcula
       el rango [0, n_local) sobre él. */
    if (n_local + n_ghost > *capacity) {
        *capacity = n_local + n_ghost;
        *p = realloc(*p, (size_t)(*capacity) * sizeof(Particle));
    }

    MPI_Alltoallv(send.buf, scount, sdispl, ptype,
                  *p + n_local, rcount, rdispl, ptype, comm);

    free(send.buf);
    free(boxes);
    free(scount); free(sdispl); free(rcount); free(rdispl);

    return n_ghost;
}

double let_mass_error(const Particle *p, int n_local, int n_ghost,
                      MPI_Comm comm) {
    double m_local = 0.0, m_seen = 0.0;
    for (int i = 0; i < n_local; i++)           m_local += p[i].mass;
    for (int i = 0; i < n_local + n_ghost; i++) m_seen  += p[i].mass;

    double m_total;
    MPI_Allreduce(&m_local, &m_total, 1, MPI_DOUBLE, MPI_SUM, comm);

    double err = (m_total > 0.0) ? fabs(m_seen - m_total) / m_total : 0.0;
    double err_max;
    MPI_Allreduce(&err, &err_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    return err_max;
}

#endif /* USE_MPI */
