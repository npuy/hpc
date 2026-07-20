#include "leapfrog.h"
#include "vec3.h"

/* Ambos operadores son elemento a elemento: sin dependencias entre iteraciones
   y con costo uniforme, por lo que schedule(static) es óptimo. */
void leapfrog_kick(Particle *p, int n, double dt_half) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        p[i].vel[0] += p[i].acc[0] * dt_half;
        p[i].vel[1] += p[i].acc[1] * dt_half;
        p[i].vel[2] += p[i].acc[2] * dt_half;
    }
}

void leapfrog_drift(Particle *p, int n, double dt) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        p[i].pos[0] += p[i].vel[0] * dt;
        p[i].pos[1] += p[i].vel[1] * dt;
        p[i].pos[2] += p[i].vel[2] * dt;
    }
}
