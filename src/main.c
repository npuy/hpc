#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "force.h"
#include "leapfrog.h"
#include "morton.h"
#include "validation.h"
#include "vec3.h"

static void print_usage(const char *prog) {
    printf("Uso: %s [opciones]\n", prog);
    printf("  -n N          Numero de particulas (default: 1000)\n");
    printf("  -dt DT        Paso temporal (default: 0.001)\n");
    printf("  -t T_END      Tiempo final (default: 1.0)\n");
    printf("  -s EPS        Softening (default: 0.01)\n");
    printf("  -o FILE       Archivo de salida CSV\n");
    printf("  --init TYPE   Tipo de init: plummer|uniform|twobody (default: plummer)\n");
    printf("  --validate    Correr suite de validacion\n");
    printf("  --seed SEED   Semilla aleatoria (default: 42)\n");
}

int main(int argc, char *argv[]) {
    int N = 1000;
    double dt = 0.001;
    double t_end = 1.0;
    double softening = 0.01;
    char *output = NULL;
    char *init_type = "plummer";
    int validate = 0;
    unsigned seed = 42;
    int snapshot_interval = 100;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i+1 < argc)
            N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-dt") == 0 && i+1 < argc)
            dt = atof(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i+1 < argc)
            t_end = atof(argv[++i]);
        else if (strcmp(argv[i], "-s") == 0 && i+1 < argc)
            softening = atof(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc)
            output = argv[++i];
        else if (strcmp(argv[i], "--init") == 0 && i+1 < argc)
            init_type = argv[++i];
        else if (strcmp(argv[i], "--validate") == 0)
            validate = 1;
        else if (strcmp(argv[i], "--seed") == 0 && i+1 < argc)
            seed = (unsigned)atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }

    if (validate)
        return run_all_tests();

    printf("N-Body Simulacion Secuencial\n");
    printf("  N = %d, dt = %g, t_end = %g, softening = %g\n", N, dt, t_end, softening);
    printf("  Init: %s, seed: %u\n\n", init_type, seed);

    Particle *particles = malloc(N * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Error: no se pudo asignar memoria para %d particulas\n", N);
        return 1;
    }

    if (strcmp(init_type, "twobody") == 0) {
        init_two_body(particles, &N);
    } else if (strcmp(init_type, "uniform") == 0) {
        init_uniform_cube(particles, N, seed);
    } else {
        init_plummer(particles, N, seed);
    }

    /* Morton sort initial */
    double mn[3] = {1e30, 1e30, 1e30};
    double mx[3] = {-1e30, -1e30, -1e30};
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            if (particles[i].pos[k] < mn[k]) mn[k] = particles[i].pos[k];
            if (particles[i].pos[k] > mx[k]) mx[k] = particles[i].pos[k];
        }
    }
    morton_sort(particles, N, mn, mx);

    compute_forces_direct(particles, N, softening);

    double E0 = kinetic_energy(particles, N) + potential_energy(particles, N, softening);
    printf("Energia inicial E0 = %.10e\n", E0);

    if (output)
        write_particles_csv(output, particles, N, 0.0);

    struct timespec ts_start, ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int step = 0;
    int total_steps = (int)(t_end / dt);

    for (double t = 0; t < t_end - dt * 0.5; t += dt, step++) {
        leapfrog_kick(particles, N, dt / 2.0);
        leapfrog_drift(particles, N, dt);
        compute_forces_direct(particles, N, softening);
        leapfrog_kick(particles, N, dt / 2.0);

        if (snapshot_interval > 0 && step % snapshot_interval == 0) {
            double E = kinetic_energy(particles, N) + potential_energy(particles, N, softening);
            double drift = fabs(E - E0) / fabs(E0);
            printf("  paso %d/%d  t=%.4f  E=%.10e  drift=%.4e\n",
                   step, total_steps, t + dt, E, drift);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double wall_time = (ts_end.tv_sec - ts_start.tv_sec)
                     + (ts_end.tv_nsec - ts_start.tv_nsec) * 1e-9;

    double Ef = kinetic_energy(particles, N) + potential_energy(particles, N, softening);
    double final_drift = fabs(Ef - E0) / fabs(E0);

    printf("\n=== Resultado final ===\n");
    printf("  Pasos: %d\n", step);
    printf("  Energia final Ef = %.10e\n", Ef);
    printf("  Energy drift: %.6e\n", final_drift);
    printf("  Tiempo de ejecucion: %.3f s\n", wall_time);
    printf("  Tiempo por paso: %.3f ms\n", (wall_time / step) * 1000.0);

    if (output) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s_final.csv", output);
        write_particles_csv(fname, particles, N, t_end);
        printf("  Snapshot final: %s\n", fname);
    }

    free(particles);
    return 0;
}
