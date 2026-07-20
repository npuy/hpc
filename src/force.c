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

/*
 * Variante paralelizable: NO usa la tercera ley de Newton. Cada iteracion i
 * recorre todos los j != i y escribe unicamente en p[i].acc, por lo que no hay
 * carrera de datos y el resultado es bit a bit reproducible con cualquier
 * numero de hilos.
 *
 * compute_forces_direct acumula en p[j].acc dentro del bucle interno; con
 * "parallel for" sobre i, dos hilos escribirian el mismo p[j].acc. Esa version
 * queda intacta como oraculo secuencial de los tests 1-8.
 *
 * Costo: 2x los flops de la version simetrica (N(N-1) pares en vez de N(N-1)/2),
 * compensado con creces por el paralelismo a partir de 2 hilos.
 */
void compute_forces_direct_par(Particle *p, int n, double softening) {
    double eps2 = softening * softening;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double acc[3] = {0, 0, 0};
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double r[3];
            VEC3_SUB(r, p[j].pos, p[i].pos);

            double dist2 = VEC3_DOT(r, r) + eps2;
            double inv_dist = 1.0 / sqrt(dist2);
            double inv_dist3 = inv_dist * inv_dist * inv_dist;

            double f = p[j].mass * inv_dist3;
            acc[0] += f * r[0];
            acc[1] += f * r[1];
            acc[2] += f * r[2];
        }
        VEC3_COPY(p[i].acc, acc);
    }
}

double kinetic_energy(const Particle *p, int n) {
    return kinetic_energy_range(p, n, 0, n);
}

double kinetic_energy_range(const Particle *p, int n, int i0, int i1) {
    (void)n;
    double ke = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ke)
    for (int i = i0; i < i1; i++)
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
