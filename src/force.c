#include "force.h"
#include "vec3.h"
#include <math.h>
#include <string.h>

void compute_forces_direct(Particle *p, int n, double softening) {
    double eps2 = softening * softening;

    for (int i = 0; i < n; i++)
        VEC3_ZERO(p[i].acc);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double r[3];
            VEC3_SUB(r, p[j].pos, p[i].pos);

            double dist2 = VEC3_DOT(r, r) + eps2;
            double inv_dist = 1.0 / sqrt(dist2);
            double inv_dist3 = inv_dist * inv_dist * inv_dist;

            double fi = p[j].mass * inv_dist3;
            double fj = p[i].mass * inv_dist3;

            p[i].acc[0] += fi * r[0];
            p[i].acc[1] += fi * r[1];
            p[i].acc[2] += fi * r[2];

            p[j].acc[0] -= fj * r[0];
            p[j].acc[1] -= fj * r[1];
            p[j].acc[2] -= fj * r[2];
        }
    }
}

double kinetic_energy(const Particle *p, int n) {
    double ke = 0.0;
    for (int i = 0; i < n; i++)
        ke += 0.5 * p[i].mass * VEC3_DOT(p[i].vel, p[i].vel);
    return ke;
}

double potential_energy(const Particle *p, int n, double softening) {
    double pe = 0.0;
    double eps2 = softening * softening;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double r[3];
            VEC3_SUB(r, p[j].pos, p[i].pos);
            double dist = sqrt(VEC3_DOT(r, r) + eps2);
            pe -= p[i].mass * p[j].mass / dist;
        }
    }
    return pe;
}
