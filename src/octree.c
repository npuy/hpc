#include "octree.h"
#include "vec3.h"
#include <stdlib.h>

#define MAX_DEPTH 21       /* límite de resolución (coincide con Morton) */
#define INIT_CAPACITY 1024 /* nodos iniciales del pool */

/*
 * Reserva un nodo nuevo en el pool (creciendo por duplicación) y lo inicializa
 * como hoja vacía. Retorna su índice. Puede reubicar t->nodes vía realloc.
 */
static int alloc_node(Octree *t) {
    if (t->count == t->capacity) {
        t->capacity *= 2;
        t->nodes = realloc(t->nodes, t->capacity * sizeof(OctreeNode));
    }
    int idx = t->count++;
    OctreeNode *nd = &t->nodes[idx];
    nd->mass = 0.0;
    VEC3_ZERO(nd->cm);
    for (int k = 0; k < 8; k++) nd->children[k] = -1;
    nd->particle = -1;
    nd->is_leaf = 1;
    return idx;
}

void octree_bounds(const Particle *p, int n, double min[3], double max[3]) {
    double mn[3] = { 1e30,  1e30,  1e30};
    double mx[3] = {-1e30, -1e30, -1e30};
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            if (p[i].pos[k] < mn[k]) mn[k] = p[i].pos[k];
            if (p[i].pos[k] > mx[k]) mx[k] = p[i].pos[k];
        }
    }

    /* lado = mayor extensión de los 3 ejes -> cubo */
    double side = 0.0;
    for (int k = 0; k < 3; k++)
        if (mx[k] - mn[k] > side) side = mx[k] - mn[k];
    if (side <= 0.0) side = 1.0;
    side += 2.0 * (side * 0.01);  /* margen 1% por lado */

    for (int k = 0; k < 3; k++) {
        double c = 0.5 * (mn[k] + mx[k]);
        min[k] = c - 0.5 * side;
        max[k] = c + 0.5 * side;
    }
}

/* Octante (0..7) del punto pos respecto al centro de la celda del nodo idx. */
static int octant_of(const Octree *t, int idx, const double pos[3]) {
    const OctreeNode *nd = &t->nodes[idx];
    double cx = 0.5 * (nd->min[0] + nd->max[0]);
    double cy = 0.5 * (nd->min[1] + nd->max[1]);
    double cz = 0.5 * (nd->min[2] + nd->max[2]);
    return (pos[0] >= cx ? 1 : 0) |
           (pos[1] >= cy ? 2 : 0) |
           (pos[2] >= cz ? 4 : 0);
}

/*
 * Convierte la hoja idx en nodo interno: crea sus 8 hijos como subcubos y la
 * marca is_leaf=0. Captura el bbox por valor porque alloc_node puede reubicar
 * el pool.
 */
static void subdivide(Octree *t, int idx) {
    double mn[3], mx[3];
    VEC3_COPY(mn, t->nodes[idx].min);
    VEC3_COPY(mx, t->nodes[idx].max);
    double c[3] = {0.5*(mn[0]+mx[0]), 0.5*(mn[1]+mx[1]), 0.5*(mn[2]+mx[2])};

    for (int ch = 0; ch < 8; ch++) {
        int cidx = alloc_node(t);            /* puede realloc: refetch abajo */
        OctreeNode *nd = &t->nodes[cidx];
        nd->min[0] = (ch & 1) ? c[0] : mn[0];
        nd->max[0] = (ch & 1) ? mx[0] : c[0];
        nd->min[1] = (ch & 2) ? c[1] : mn[1];
        nd->max[1] = (ch & 2) ? mx[1] : c[1];
        nd->min[2] = (ch & 4) ? c[2] : mn[2];
        nd->max[2] = (ch & 4) ? mx[2] : c[2];
        t->nodes[idx].children[ch] = cidx;
    }
    t->nodes[idx].is_leaf = 0;
}

/* Inserta la partícula i en el subárbol enraizado en idx (índice fresco). */
static void octree_insert(Octree *t, int idx, const Particle *p, int i, int depth) {
    if (t->nodes[idx].is_leaf) {
        if (t->nodes[idx].particle == -1) {   /* hoja vacía */
            t->nodes[idx].particle = i;
            return;
        }
        if (depth >= MAX_DEPTH) {             /* coincidentes: encadenar */
            t->next[i] = t->nodes[idx].particle;
            t->nodes[idx].particle = i;
            return;
        }
        /* hoja ocupada: subdividir y reubicar la partícula existente j */
        int j = t->nodes[idx].particle;
        t->nodes[idx].particle = -1;
        subdivide(t, idx);
        int oj = octant_of(t, idx, p[j].pos);
        octree_insert(t, t->nodes[idx].children[oj], p, j, depth + 1);
    }
    /* nodo interno: descender al octante de i */
    int oi = octant_of(t, idx, p[i].pos);
    octree_insert(t, t->nodes[idx].children[oi], p, i, depth + 1);
}

Octree *octree_build_box(const Particle *p, int n,
                         const double min[3], const double max[3]) {
    Octree *t = malloc(sizeof(Octree));
    t->capacity = INIT_CAPACITY;
    t->nodes = malloc(t->capacity * sizeof(OctreeNode));
    t->count = 0;
    t->next = malloc((n > 0 ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) t->next[i] = -1;

    t->root = alloc_node(t);
    VEC3_COPY(t->nodes[t->root].min, min);
    VEC3_COPY(t->nodes[t->root].max, max);
    t->size = max[0] - min[0];

    for (int i = 0; i < n; i++)
        octree_insert(t, t->root, p, i, 0);

    return t;
}

Octree *octree_build(const Particle *p, int n) {
    double mn[3], mx[3];
    octree_bounds(p, n, mn, mx);
    return octree_build_box(p, n, mn, mx);
}

void octree_insert_particles(Octree *t, const Particle *p, int n_old, int n_new) {
    if (n_new <= n_old) return;

    /* next[] está indexado por índice de partícula: hay que agrandarlo. realloc
       preserva las cadenas ya armadas para las primeras n_old. */
    t->next = realloc(t->next, n_new * sizeof(int));
    for (int i = n_old; i < n_new; i++) t->next[i] = -1;

    for (int i = n_old; i < n_new; i++)
        octree_insert(t, t->root, p, i, 0);
}

/* Post-orden: acumula masa y centro de masa desde las hojas hacia la raíz. */
static void compute_mass_rec(Octree *t, int idx, const Particle *p) {
    if (t->nodes[idx].is_leaf) {
        double m = 0.0, cm[3] = {0, 0, 0};
        for (int k = t->nodes[idx].particle; k != -1; k = t->next[k]) {
            m += p[k].mass;
            cm[0] += p[k].mass * p[k].pos[0];
            cm[1] += p[k].mass * p[k].pos[1];
            cm[2] += p[k].mass * p[k].pos[2];
        }
        t->nodes[idx].mass = m;
        if (m > 0.0) {
            t->nodes[idx].cm[0] = cm[0] / m;
            t->nodes[idx].cm[1] = cm[1] / m;
            t->nodes[idx].cm[2] = cm[2] / m;
        }
        return;
    }

    double m = 0.0, cm[3] = {0, 0, 0};
    for (int ch = 0; ch < 8; ch++) {
        int c = t->nodes[idx].children[ch];
        if (c == -1) continue;
        compute_mass_rec(t, c, p);
        double cmass = t->nodes[c].mass;
        m += cmass;
        cm[0] += cmass * t->nodes[c].cm[0];
        cm[1] += cmass * t->nodes[c].cm[1];
        cm[2] += cmass * t->nodes[c].cm[2];
    }
    t->nodes[idx].mass = m;
    if (m > 0.0) {
        t->nodes[idx].cm[0] = cm[0] / m;
        t->nodes[idx].cm[1] = cm[1] / m;
        t->nodes[idx].cm[2] = cm[2] / m;
    }
}

void octree_compute_mass(Octree *t, const Particle *p) {
    compute_mass_rec(t, t->root, p);
}

void octree_free(Octree *t) {
    if (!t) return;
    free(t->nodes);
    free(t->next);
    free(t);
}
