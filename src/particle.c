#include "particle.h"
#include "vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void init_two_body(Particle *p, int *n) {
    *n = 2;
    memset(p, 0, 2 * sizeof(Particle));

    double M = 1.0;
    double m = M / 2.0;
    double sep = 1.0;
    double r = sep / 2.0;
    /* v for circular orbit: v = sqrt(G*m / (4*r)) with G=1, total_mass=1 */
    double v = sqrt(M / (4.0 * sep));

    p[0].mass = m;
    p[0].pos[0] = -r;
    p[0].vel[1] = -v;

    p[1].mass = m;
    p[1].pos[0] = r;
    p[1].vel[1] = v;
}

void init_plummer(Particle *p, int n, unsigned seed) {
    srand(seed);
    double mass_each = 1.0 / n;
    double a = 1.0; /* Plummer scale radius */

    for (int i = 0; i < n; i++) {
        memset(&p[i], 0, sizeof(Particle));
        p[i].mass = mass_each;

        /* Plummer model: radius from inverse CDF */
        double X1 = ((double)rand() / RAND_MAX);
        if (X1 < 1e-10) X1 = 1e-10;
        double ri = a / sqrt(pow(X1, -2.0/3.0) - 1.0);

        /* Uniform direction on sphere */
        double costheta = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        double sintheta = sqrt(1.0 - costheta * costheta);
        double phi = 2.0 * M_PI * ((double)rand() / RAND_MAX);

        p[i].pos[0] = ri * sintheta * cos(phi);
        p[i].pos[1] = ri * sintheta * sin(phi);
        p[i].pos[2] = ri * costheta;

        /* Velocity: rejection sampling (Aarseth, Henon, Wielen 1974) */
        double ve = sqrt(2.0) * pow(1.0 + ri * ri / (a * a), -0.25);
        double q, g;
        do {
            q = ((double)rand() / RAND_MAX);
            g = q * q * pow(1.0 - q * q, 3.5);
        } while (10.0 * ((double)rand() / RAND_MAX) > g / 0.1);

        double vi = q * ve;

        costheta = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        sintheta = sqrt(1.0 - costheta * costheta);
        phi = 2.0 * M_PI * ((double)rand() / RAND_MAX);

        p[i].vel[0] = vi * sintheta * cos(phi);
        p[i].vel[1] = vi * sintheta * sin(phi);
        p[i].vel[2] = vi * costheta;
    }

    /* Center of mass correction */
    double cm_pos[3] = {0, 0, 0};
    double cm_vel[3] = {0, 0, 0};
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            cm_pos[k] += p[i].mass * p[i].pos[k];
            cm_vel[k] += p[i].mass * p[i].vel[k];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            p[i].pos[k] -= cm_pos[k];
            p[i].vel[k] -= cm_vel[k];
        }
    }
}

void init_uniform_cube(Particle *p, int n, unsigned seed) {
    srand(seed);
    double mass_each = 1.0 / n;

    for (int i = 0; i < n; i++) {
        memset(&p[i], 0, sizeof(Particle));
        p[i].mass = mass_each;
        p[i].pos[0] = (double)rand() / RAND_MAX - 0.5;
        p[i].pos[1] = (double)rand() / RAND_MAX - 0.5;
        p[i].pos[2] = (double)rand() / RAND_MAX - 0.5;
    }
}

void write_particles_csv(const char *fname, const Particle *p, int n, double t) {
    FILE *f = fopen(fname, "w");
    if (!f) {
        fprintf(stderr, "Error: no se pudo abrir %s para escritura\n", fname);
        return;
    }
    fprintf(f, "t,mass,x,y,z,vx,vy,vz\n");
    for (int i = 0; i < n; i++) {
        fprintf(f, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                t, p[i].mass,
                p[i].pos[0], p[i].pos[1], p[i].pos[2],
                p[i].vel[0], p[i].vel[1], p[i].vel[2]);
    }
    fclose(f);
}

int read_particles_csv(const char *fname, Particle *p, int *n) {
    FILE *f = fopen(fname, "r");
    if (!f) return -1;

    char line[512];
    /* skip header */
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    int count = 0;
    double t;
    while (fgets(line, sizeof(line), f)) {
        Particle *pi = &p[count];
        memset(pi, 0, sizeof(Particle));
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &t, &pi->mass,
                   &pi->pos[0], &pi->pos[1], &pi->pos[2],
                   &pi->vel[0], &pi->vel[1], &pi->vel[2]) == 8) {
            count++;
        }
    }
    *n = count;
    fclose(f);
    return 0;
}
