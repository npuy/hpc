#include "barnes_hut.h"
#include "vec3.h"
#include <math.h>

/*
 * Acumula en acc[] la aceleración que el subárbol idx ejerce sobre la
 * partícula i. eps2 = softening². theta2 = theta² (para comparar sin raíces).
 * El árbol es inmutable durante el cálculo, por lo que cachear el puntero al
 * nodo es seguro (no hay realloc).
 */
static void accumulate_force(const Octree *t, int idx, int i, const Particle *p,
                             double theta2, double eps2, double acc[3]) {
    const OctreeNode *nd = &t->nodes[idx];
    if (nd->mass == 0.0) return;   /* celda vacía */

    if (nd->is_leaf) {
        /* interacción directa con cada partícula de la hoja, salteando i */
        for (int k = nd->particle; k != -1; k = t->next[k]) {
            if (k == i) continue;
            double r[3];
            VEC3_SUB(r, p[k].pos, p[i].pos);
            double dist2 = VEC3_DOT(r, r) + eps2;
            double inv = 1.0 / sqrt(dist2);
            double inv3 = inv * inv * inv;
            double f = p[k].mass * inv3;
            acc[0] += f * r[0];
            acc[1] += f * r[1];
            acc[2] += f * r[2];
        }
        return;
    }

    /* nodo interno: criterio de apertura s/d < theta  <=>  s² < theta²·d² */
    double r[3];
    VEC3_SUB(r, nd->cm, p[i].pos);
    double dist2 = VEC3_DOT(r, r) + eps2;
    double s = nd->max[0] - nd->min[0];

    if (s * s < theta2 * dist2) {
        double inv = 1.0 / sqrt(dist2);
        double inv3 = inv * inv * inv;
        double f = nd->mass * inv3;
        acc[0] += f * r[0];
        acc[1] += f * r[1];
        acc[2] += f * r[2];
    } else {
        for (int ch = 0; ch < 8; ch++) {
            int c = nd->children[ch];
            if (c != -1)
                accumulate_force(t, c, i, p, theta2, eps2, acc);
        }
    }
}

void compute_forces_bh(const Octree *t, Particle *p, int n,
                       double theta, double softening) {
    double eps2 = softening * softening;
    double theta2 = theta * theta;

    for (int i = 0; i < n; i++) {
        double acc[3] = {0, 0, 0};
        accumulate_force(t, t->root, i, p, theta2, eps2, acc);
        VEC3_COPY(p[i].acc, acc);
    }
}

/* Potencial en la partícula i por recorrido del subárbol idx. */
static double potential_at(const Octree *t, int idx, int i, const Particle *p,
                           double theta2, double eps2) {
    const OctreeNode *nd = &t->nodes[idx];
    if (nd->mass == 0.0) return 0.0;

    if (nd->is_leaf) {
        double phi = 0.0;
        for (int k = nd->particle; k != -1; k = t->next[k]) {
            if (k == i) continue;
            double r[3];
            VEC3_SUB(r, p[k].pos, p[i].pos);
            double dist = sqrt(VEC3_DOT(r, r) + eps2);
            phi -= p[k].mass / dist;
        }
        return phi;
    }

    double r[3];
    VEC3_SUB(r, nd->cm, p[i].pos);
    double dist2 = VEC3_DOT(r, r) + eps2;
    double s = nd->max[0] - nd->min[0];

    if (s * s < theta2 * dist2)
        return -nd->mass / sqrt(dist2);

    double phi = 0.0;
    for (int ch = 0; ch < 8; ch++) {
        int c = nd->children[ch];
        if (c != -1)
            phi += potential_at(t, c, i, p, theta2, eps2);
    }
    return phi;
}

double potential_energy_bh(const Octree *t, const Particle *p, int n,
                           double theta, double softening) {
    double eps2 = softening * softening;
    double theta2 = theta * theta;
    double U = 0.0;
    for (int i = 0; i < n; i++)
        U += 0.5 * p[i].mass * potential_at(t, t->root, i, p, theta2, eps2);
    return U;
}
