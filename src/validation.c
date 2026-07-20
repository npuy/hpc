#include "validation.h"
#include "particle.h"
#include "force.h"
#include "leapfrog.h"
#include "morton.h"
#include "octree.h"
#include "barnes_hut.h"
#include "vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "mpi_types.h"
#include "domain.h"
#include "migration.h"
#endif

int test_two_body_orbit(void) {
    printf("=== Test 1: Orbita circular de dos cuerpos ===\n");

    Particle p[2];
    int n;
    init_two_body(p, &n);

    double pos0_x = p[0].pos[0];
    double pos0_y = p[0].pos[1];
    double pos1_x = p[1].pos[0];
    double pos1_y = p[1].pos[1];

    double softening = 0.001;
    double sep = p[1].pos[0] - p[0].pos[0];
    /* Period = 2*pi*r^(3/2) / sqrt(G*M) with r=sep, G=1, M=1 */
    double T = 2.0 * M_PI * pow(sep, 1.5);
    double dt = T / 10000.0;
    int steps = (int)(T / dt);

    compute_forces_direct(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double err0 = sqrt((p[0].pos[0]-pos0_x)*(p[0].pos[0]-pos0_x) +
                        (p[0].pos[1]-pos0_y)*(p[0].pos[1]-pos0_y));
    double err1 = sqrt((p[1].pos[0]-pos1_x)*(p[1].pos[0]-pos1_x) +
                        (p[1].pos[1]-pos1_y)*(p[1].pos[1]-pos1_y));
    double max_err = fmax(err0, err1);
    double ref = sep;

    printf("  Periodo orbital T = %.4f, dt = %.6f, pasos = %d\n", T, dt, steps);
    printf("  Error pos cuerpo 0: %.6e\n", err0);
    printf("  Error pos cuerpo 1: %.6e\n", err1);
    printf("  Error relativo max: %.4f%%\n", (max_err / ref) * 100.0);

    int pass = (max_err / ref) < 0.01;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    return pass;
}

int test_energy_conservation(void) {
    printf("=== Test 2: Conservacion de energia (N=100, 100 pasos) ===\n");

    int n = 100;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double softening = 0.01;
    double dt = 0.001;
    int steps = 100;

    compute_forces_direct(p, n, softening);
    double E0 = kinetic_energy(p, n) + potential_energy(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double Ef = kinetic_energy(p, n) + potential_energy(p, n, softening);
    double drift = fabs(Ef - E0) / fabs(E0);

    printf("  E0 = %.10e\n", E0);
    printf("  Ef = %.10e\n", Ef);
    printf("  Drift relativo: %.6e\n", drift);

    int pass = drift < 1e-3;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(p);
    return pass;
}

int test_momentum_conservation(void) {
    printf("=== Test 3: Conservacion de momento lineal ===\n");

    int n = 50;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 123);

    double softening = 0.01;
    double dt = 0.001;
    int steps = 50;

    double p0[3] = {0, 0, 0};
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            p0[k] += p[i].mass * p[i].vel[k];

    compute_forces_direct(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double pf[3] = {0, 0, 0};
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            pf[k] += p[i].mass * p[i].vel[k];

    double dp = sqrt((pf[0]-p0[0])*(pf[0]-p0[0]) +
                      (pf[1]-p0[1])*(pf[1]-p0[1]) +
                      (pf[2]-p0[2])*(pf[2]-p0[2]));

    printf("  Momento inicial: (%.6e, %.6e, %.6e)\n", p0[0], p0[1], p0[2]);
    printf("  Momento final:   (%.6e, %.6e, %.6e)\n", pf[0], pf[1], pf[2]);
    printf("  |dp|: %.6e\n", dp);

    int pass = dp < 1e-12;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(p);
    return pass;
}

int test_morton_ordering(void) {
    printf("=== Test 4: Morton ordering (8 esquinas del cubo) ===\n");

    Particle p[8];
    memset(p, 0, sizeof(p));
    double corners[8][3] = {
        {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
        {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}
    };
    for (int i = 0; i < 8; i++) {
        p[i].pos[0] = corners[i][0];
        p[i].pos[1] = corners[i][1];
        p[i].pos[2] = corners[i][2];
        p[i].mass = 1.0;
    }

    double mn[3] = {0, 0, 0};
    double mx[3] = {1, 1, 1};

    for (int i = 0; i < 8; i++) {
        p[i].morton = morton_encode(p[i].pos[0], p[i].pos[1], p[i].pos[2],
                                    mn, mx, 21);
    }

    printf("  Claves Morton antes de sort:\n");
    for (int i = 0; i < 8; i++)
        printf("    (%.0f,%.0f,%.0f) -> 0x%016llx\n",
               p[i].pos[0], p[i].pos[1], p[i].pos[2],
               (unsigned long long)p[i].morton);

    morton_sort(p, 8, mn, mx);

    printf("  Claves Morton despues de sort:\n");
    for (int i = 0; i < 8; i++)
        printf("    (%.0f,%.0f,%.0f) -> 0x%016llx\n",
               p[i].pos[0], p[i].pos[1], p[i].pos[2],
               (unsigned long long)p[i].morton);

    int pass = 1;
    for (int i = 1; i < 8; i++) {
        if (p[i].morton < p[i-1].morton) {
            pass = 0;
            break;
        }
    }

    /* Verify all 8 distinct keys */
    for (int i = 0; i < 8 && pass; i++)
        for (int j = i+1; j < 8 && pass; j++)
            if (p[i].morton == p[j].morton)
                pass = 0;

    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    return pass;
}

/* Error relativo L2 global entre dos campos de aceleracion. */
static double accel_l2_error(const Particle *ref, const Particle *bh, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        double d[3];
        VEC3_SUB(d, bh[i].acc, ref[i].acc);
        num += VEC3_DOT(d, d);
        den += VEC3_DOT(ref[i].acc, ref[i].acc);
    }
    return sqrt(num / den);
}

int test_tree_mass(void) {
    printf("=== Test 5: Conservacion de masa del arbol ===\n");

    int n = 1000;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double total = 0.0;
    for (int i = 0; i < n; i++) total += p[i].mass;

    Octree *t = octree_build(p, n);
    octree_compute_mass(t, p);
    double root_mass = t->nodes[t->root].mass;

    double err = fabs(root_mass - total) / total;
    printf("  Masa total (suma directa): %.15e\n", total);
    printf("  Masa de la raiz:           %.15e\n", root_mass);
    printf("  Error relativo: %.6e\n", err);

    int pass = err < 1e-12;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(p);
    return pass;
}

int test_tree_cm(void) {
    printf("=== Test 6: Centro de masa del arbol ===\n");

    int n = 1000;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double total = 0.0, cm[3] = {0, 0, 0};
    for (int i = 0; i < n; i++) {
        total += p[i].mass;
        for (int k = 0; k < 3; k++) cm[k] += p[i].mass * p[i].pos[k];
    }
    for (int k = 0; k < 3; k++) cm[k] /= total;

    Octree *t = octree_build(p, n);
    octree_compute_mass(t, p);
    double *rcm = t->nodes[t->root].cm;

    double d[3];
    VEC3_SUB(d, rcm, cm);
    double err = VEC3_NORM(d);
    printf("  CM directo: (%.6e, %.6e, %.6e)\n", cm[0], cm[1], cm[2]);
    printf("  CM raiz:    (%.6e, %.6e, %.6e)\n", rcm[0], rcm[1], rcm[2]);
    printf("  |diferencia|: %.6e\n", err);

    int pass = err < 1e-10;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(p);
    return pass;
}

int test_bh_force_error(void) {
    printf("=== Test 7: Error de fuerza Barnes-Hut vs directo (theta=0.5) ===\n");

    int n = 1000;
    double softening = 0.01, theta = 0.5;
    Particle *ref = malloc(n * sizeof(Particle));
    init_plummer(ref, n, 42);

    Particle *bh = malloc(n * sizeof(Particle));
    memcpy(bh, ref, n * sizeof(Particle));

    compute_forces_direct(ref, n, softening);

    Octree *t = octree_build(bh, n);
    octree_compute_mass(t, bh);
    compute_forces_bh(t, bh, n, theta, softening);

    double err = accel_l2_error(ref, bh, n);
    printf("  Error relativo L2 de la aceleracion: %.6e\n", err);

    int pass = err < 1e-2;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(bh);
    free(ref);
    return pass;
}

int test_bh_theta_convergence(void) {
    printf("=== Test 8: Convergencia en theta ===\n");

    int n = 1000;
    double softening = 0.01;
    double thetas[4] = {0.8, 0.5, 0.2, 0.0};

    Particle *ref = malloc(n * sizeof(Particle));
    init_plummer(ref, n, 42);
    compute_forces_direct(ref, n, softening);

    double errs[4];
    for (int idx = 0; idx < 4; idx++) {
        Particle *bh = malloc(n * sizeof(Particle));
        memcpy(bh, ref, n * sizeof(Particle));
        Octree *t = octree_build(bh, n);
        octree_compute_mass(t, bh);
        compute_forces_bh(t, bh, n, thetas[idx], softening);
        errs[idx] = accel_l2_error(ref, bh, n);
        printf("  theta = %.1f  ->  error L2 = %.6e\n", thetas[idx], errs[idx]);
        octree_free(t);
        free(bh);
    }

    int monotone = 1;
    for (int idx = 1; idx < 4; idx++)
        if (errs[idx] > errs[idx-1]) monotone = 0;

    int exact = errs[3] < 1e-9;  /* theta=0 reproduce el metodo directo */

    printf("  Monotono decreciente: %s\n", monotone ? "si" : "no");
    printf("  theta=0 reproduce directo (err < 1e-9): %s\n", exact ? "si" : "no");

    int pass = monotone && exact;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(ref);
    return pass;
}

int test_openmp_determinism(void) {
    printf("=== Test 9: Determinismo de OpenMP ===\n");

#ifndef _OPENMP
    printf("  Compilado sin OpenMP: nada que verificar.\n");
    printf("  Resultado: SKIP\n\n");
    return 1;
#else
    int n = 5000;
    double softening = 0.01, theta = 0.5;
    int thread_counts[4] = {1, 2, 4, 8};

    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);
    double mn[3], mx[3];
    octree_bounds(p, n, mn, mx);
    morton_sort(p, n, mn, mx);

    Octree *t = octree_build(p, n);
    octree_compute_mass(t, p);

    /* El arbol es inmutable durante el calculo y cada iteracion escribe solo
       p[i].acc, asi que el resultado debe ser bit a bit identico. Cualquier
       diferencia delata una carrera de datos. */
    double *ref_bh  = malloc(3 * n * sizeof(double));
    double *ref_dir = malloc(3 * n * sizeof(double));
    int pass_bh = 1, pass_dir = 1;

    for (int k = 0; k < 4; k++) {
        omp_set_num_threads(thread_counts[k]);

        compute_forces_bh(t, p, n, theta, softening);
        if (k == 0) {
            for (int i = 0; i < n; i++) VEC3_COPY(&ref_bh[3*i], p[i].acc);
        } else {
            for (int i = 0; i < n; i++)
                if (memcmp(&ref_bh[3*i], p[i].acc, 3 * sizeof(double)) != 0)
                    pass_bh = 0;
        }

        compute_forces_direct_par(p, n, softening);
        if (k == 0) {
            for (int i = 0; i < n; i++) VEC3_COPY(&ref_dir[3*i], p[i].acc);
        } else {
            for (int i = 0; i < n; i++)
                if (memcmp(&ref_dir[3*i], p[i].acc, 3 * sizeof(double)) != 0)
                    pass_dir = 0;
        }

        printf("  %d hilo(s): bh %s, directo-par %s\n", thread_counts[k],
               (k == 0) ? "(referencia)" : (pass_bh  ? "identico" : "DIFIERE"),
               (k == 0) ? "(referencia)" : (pass_dir ? "identico" : "DIFIERE"));
    }

    omp_set_num_threads(omp_get_max_threads());

    int pass = pass_bh && pass_dir;
    printf("  Igualdad bit a bit con 1/2/4/8 hilos: %s\n", pass ? "si" : "no");
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(ref_bh); free(ref_dir);
    octree_free(t);
    free(p);
    return pass;
#endif
}

#ifdef USE_MPI

/*
 * Monta el estado distribuido inicial: el rango 0 genera N particulas Plummer,
 * se reparten en bloques iguales, se calcula el bounding box global, se
 * particiona por claves Morton y se migra cada particula a su dueno.
 * Retorna n_local y deja *local listo para simular.
 */
static int setup_distributed(Particle **local, int *capacity, Domain *d,
                             MPI_Datatype ptype, int N, unsigned seed,
                             MPI_Comm comm) {
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    Particle *gen = NULL;
    if (rank == 0) {
        gen = malloc(N * sizeof(Particle));
        init_plummer(gen, N, seed);
    }

    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    for (int k = 0; k < P; k++) counts[k] = N / P + (k < N % P ? 1 : 0);
    displs[0] = 0;
    for (int k = 1; k < P; k++) displs[k] = displs[k-1] + counts[k-1];

    int n_local = counts[rank];
    *capacity = (n_local > 0 ? n_local : 1);
    *local = malloc(*capacity * sizeof(Particle));
    MPI_Scatterv(gen, counts, displs, ptype, *local, n_local, ptype, 0, comm);
    free(gen); free(counts); free(displs);

    domain_global_bounds(d, *local, n_local, comm);
    double gmin[3], gmax[3];
    for (int k = 0; k < 3; k++) { gmin[k] = d->gmin[k]; gmax[k] = d->gmax[k]; }
    morton_sort(*local, n_local, gmin, gmax);
    domain_partition(d, *local, n_local, N, comm);

    return migrate_particles(local, n_local, capacity, d, ptype, comm);
}

/* Un paso completo del bucle hibrido (KDK + migracion + replica + fuerzas). */
static int distributed_step(Particle **local, int n_local, int *capacity,
                            Domain *d, MPI_Datatype ptype, Particle *all, int N,
                            double dt, double theta, double softening,
                            MPI_Comm comm) {
    leapfrog_kick(*local, n_local, dt / 2.0);
    leapfrog_drift(*local, n_local, dt);

    n_local = migrate_particles(local, n_local, capacity, d, ptype, comm);

    int offset;
    replicate_particles(*local, n_local, all, N, &offset, ptype, comm);
    Octree *t = octree_build(all, N);
    octree_compute_mass(t, all);
    compute_forces_bh_range(t, all, N, offset, offset + n_local, theta, softening);
    for (int i = 0; i < n_local; i++) VEC3_COPY((*local)[i].acc, all[offset+i].acc);
    octree_free(t);

    leapfrog_kick(*local, n_local, dt / 2.0);
    return n_local;
}

int test_mpi_particle_conservation(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 10: Conservacion de particulas (100 pasos con migracion) ===\n");

    int N = 5000, steps = 100, capacity;
    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);

    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);
    Particle *all = malloc(N * sizeof(Particle));

    int64_t worst = N;   /* peor desviacion respecto a N a lo largo de la corrida */
    for (int s = 0; s < steps; s++) {
        n_local = distributed_step(&local, n_local, &capacity, &d, ptype, all, N,
                                   0.001, 0.5, 0.01, comm);
        int64_t nl = n_local, total;
        MPI_Allreduce(&nl, &total, 1, MPI_INT64_T, MPI_SUM, comm);
        if (total != N) worst = total;
    }

    int pass = (worst == N);
    if (rank == 0) {
        printf("  N esperado: %d\n", N);
        printf("  Suma de n_local en todos los pasos: %s\n",
               pass ? "siempre N" : "DESVIACION DETECTADA");
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(all); free(local);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_mpi_identity_checksum(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 11: Identidad de particulas (sin perdidas ni duplicados) ===\n");

    int N = 5000, steps = 100, capacity;
    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);

    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);
    Particle *all = malloc(N * sizeof(Particle));

    int64_t s0, s20;
    particle_checksum(local, n_local, &s0, &s20, comm);

    for (int s = 0; s < steps; s++)
        n_local = distributed_step(&local, n_local, &capacity, &d, ptype, all, N,
                                   0.001, 0.5, 0.01, comm);

    int64_t s1, s21;
    particle_checksum(local, n_local, &s1, &s21, comm);

    int pass = (s0 == s1 && s20 == s21);
    if (rank == 0) {
        printf("  Sum(id)  inicial=%lld  final=%lld\n", (long long)s0, (long long)s1);
        printf("  Sum(id2) inicial=%lld  final=%lld\n", (long long)s20, (long long)s21);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(all); free(local);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_mpi_partition_validity(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);
    if (rank == 0)
        printf("=== Test 12: Particiones disjuntas y cubrientes ===\n");

    int N = 5000, capacity;
    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);

    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);

    /* a) splitters no decrecientes. No se exige crecimiento estricto: con claves
          repetidas dos splitters pueden coincidir y dejar un tramo vacio, que es
          un estado legal (ese proceso simplemente queda sin particulas). */
    int monotone = 1;
    for (int k = 0; k < P; k++)
        if (d.splitters[k] > d.splitters[k+1]) monotone = 0;

    /* b) cubren todo el espacio de claves */
    int covers = (d.splitters[0] == 0 && d.splitters[P] == UINT64_MAX);

    /* c) toda particula local cae dentro del tramo de su proceso */
    int inside = 1;
    for (int i = 0; i < n_local; i++)
        if (local[i].morton < d.splitters[rank] ||
            local[i].morton >= d.splitters[rank+1]) inside = 0;

    /* d) desbalance del reparto */
    int64_t nl = n_local, nmax, ntot;
    MPI_Allreduce(&nl, &nmax, 1, MPI_INT64_T, MPI_MAX, comm);
    MPI_Allreduce(&nl, &ntot, 1, MPI_INT64_T, MPI_SUM, comm);
    double imbalance = (double)nmax / ((double)ntot / P);

    int flags[3] = {monotone, covers, inside}, all_flags[3];
    MPI_Allreduce(flags, all_flags, 3, MPI_INT, MPI_MIN, comm);

    int pass = all_flags[0] && all_flags[1] && all_flags[2] && (imbalance < 1.15);
    if (rank == 0) {
        printf("  Splitters no decrecientes: %s\n", all_flags[0] ? "si" : "NO");
        printf("  Cubren [0, UINT64_MAX]:    %s\n", all_flags[1] ? "si" : "NO");
        printf("  Particulas dentro de su tramo: %s\n", all_flags[2] ? "si" : "NO");
        printf("  Desbalance max/avg: %.4f (umbral 1.15)\n", imbalance);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(local);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_mpi_forces_vs_sequential(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 13: Fuerzas distribuidas vs secuenciales ===\n");

    int N = 2000, capacity;
    double theta = 0.5, softening = 0.01;

    /* Referencia: el mismo calculo que hace la version secuencial. */
    Particle *ref = malloc(N * sizeof(Particle));
    init_plummer(ref, N, 42);
    double mn[3], mx[3];
    octree_bounds(ref, N, mn, mx);
    morton_sort(ref, N, mn, mx);
    Octree *rt = octree_build(ref, N);
    octree_compute_mass(rt, ref);
    compute_forces_bh(rt, ref, N, theta, softening);
    octree_free(rt);

    /* Camino distribuido completo: particionar, migrar, replicar, calcular el
       tramo propio. */
    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);
    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);

    Particle *all = malloc(N * sizeof(Particle));
    int offset;
    replicate_particles(local, n_local, all, N, &offset, ptype, comm);
    Octree *t = octree_build(all, N);
    octree_compute_mass(t, all);
    compute_forces_bh_range(t, all, N, offset, offset + n_local, theta, softening);
    for (int i = 0; i < n_local; i++) VEC3_COPY(local[i].acc, all[offset+i].acc);
    octree_free(t);

    /* Reunir las aceleraciones calculadas por todos los procesos. */
    replicate_particles(local, n_local, all, N, &offset, ptype, comm);

    /* La migracion permuta las particulas: el emparejamiento con la referencia
       tiene que ser por id, no por posicion en el arreglo. Para esto existe id. */
    int *pos_of_id = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) pos_of_id[ref[i].id] = i;

    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
        const Particle *a = &all[i];
        const Particle *r = &ref[pos_of_id[a->id]];
        double diff[3];
        VEC3_SUB(diff, a->acc, r->acc);
        double den = VEC3_NORM(r->acc);
        if (den > 0.0) {
            double err = VEC3_NORM(diff) / den;
            if (err > max_err) max_err = err;
        }
    }

    int pass = (max_err < 1e-12);
    if (rank == 0) {
        printf("  Error relativo maximo (emparejado por id): %.6e\n", max_err);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(pos_of_id); free(all); free(local); free(ref);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

#endif /* USE_MPI */

int run_all_tests(void) {
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int passed = 0, total = 13;

    /* Los tests 1-9 son puramente locales: los corre y los imprime el rango 0.
       Los tests 10-13 son colectivos y los deben ejecutar TODOS los procesos,
       aunque solo el rango 0 imprima. */
    if (rank == 0) {
        printf("\n========== SUITE DE VALIDACION (MPI) ==========\n\n");
        passed += test_two_body_orbit();
        passed += test_energy_conservation();
        passed += test_momentum_conservation();
        passed += test_morton_ordering();
        passed += test_tree_mass();
        passed += test_tree_cm();
        passed += test_bh_force_error();
        passed += test_bh_theta_convergence();
        passed += test_openmp_determinism();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int collective = 0;
    collective += test_mpi_particle_conservation();
    collective += test_mpi_identity_checksum();
    collective += test_mpi_partition_validity();
    collective += test_mpi_forces_vs_sequential();

    int rc = 0;
    if (rank == 0) {
        passed += collective;
        printf("========== RESULTADO: %d/%d tests pasaron ==========\n\n", passed, total);
        rc = (passed == total) ? 0 : 1;
    }
    MPI_Bcast(&rc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return rc;
#else
    printf("\n========== SUITE DE VALIDACION ==========\n\n");

    int total = 9, passed = 0;
    passed += test_two_body_orbit();
    passed += test_energy_conservation();
    passed += test_momentum_conservation();
    passed += test_morton_ordering();
    passed += test_tree_mass();
    passed += test_tree_cm();
    passed += test_bh_force_error();
    passed += test_bh_theta_convergence();
    passed += test_openmp_determinism();

    printf("========== RESULTADO: %d/%d tests pasaron ==========\n", passed, total);
    printf("(los tests 10-13 requieren MPI: mpirun -np 4 ./nbody_mpi --validate)\n\n");
    return (passed == total) ? 0 : 1;
#endif
}
