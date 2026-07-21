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
#include "let.h"
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
static int setup_distributed_ex(Particle **local, int *capacity, Domain *d,
                                MPI_Datatype ptype, int N, unsigned seed,
                                int use_orb, MPI_Comm comm) {
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
    if (use_orb) domain_partition_orb(d, *local, n_local, 0, comm);
    else         domain_partition(d, *local, n_local, N, comm);

    return migrate_particles(local, n_local, capacity, d, ptype, comm);
}

/* La descomposicion por defecto es ORB desde la semana 5. */
static int setup_distributed(Particle **local, int *capacity, Domain *d,
                             MPI_Datatype ptype, int N, unsigned seed,
                             MPI_Comm comm) {
    return setup_distributed_ex(local, capacity, d, ptype, N, seed, 1, comm);
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
        printf("=== Test 12: Descomposicion valida (ORB y Morton) ===\n");

    const int N = 5000;
    int capacity, ok_orb, ok_morton;

    /* --- a) ORB: cajas disjuntas, cubrientes, y cada particula en la suya --- */
    {
        Particle *local = NULL;
        Domain d;
        MPI_Datatype ptype;
        mpi_particle_type_init(&ptype);
        domain_init(&d, comm);
        int n_local = setup_distributed_ex(&local, &capacity, &d, ptype, N, 42, 1, comm);

        /* Disjuntas: el volumen de interseccion de cada par debe ser nulo. */
        int disjoint = 1;
        for (int a = 0; a < P && disjoint; a++) {
            double amin[3], amax[3];
            domain_box(&d, a, amin, amax);
            for (int b = a + 1; b < P && disjoint; b++) {
                double bmin[3], bmax[3];
                domain_box(&d, b, bmin, bmax);
                double v = 1.0;
                for (int k = 0; k < 3; k++) {
                    double lo = amin[k] > bmin[k] ? amin[k] : bmin[k];
                    double hi = amax[k] < bmax[k] ? amax[k] : bmax[k];
                    double e  = hi - lo;
                    v *= (e > 0.0 ? e : 0.0);
                }
                if (v > 0.0) disjoint = 0;
            }
        }

        /* Cubrientes: la suma de los volumenes reconstruye el del dominio. */
        double vsum = 0.0, vglobal = 1.0;
        for (int a = 0; a < P; a++) {
            double amin[3], amax[3];
            domain_box(&d, a, amin, amax);
            double v = 1.0;
            for (int k = 0; k < 3; k++) v *= (amax[k] - amin[k]);
            vsum += v;
        }
        for (int k = 0; k < 3; k++) vglobal *= (d.gmax[k] - d.gmin[k]);
        int covers = (fabs(vsum - vglobal) / vglobal < 1e-9);

        /* Cada particula local dentro de la caja de su proceso, y coherencia
           entre domain_owner_pos y el rank: si difieren, la migracion habria
           mandado la particula a otro lado y aca aparece. */
        double mymin[3], mymax[3];
        domain_box(&d, rank, mymin, mymax);
        int inside = 1;
        for (int i = 0; i < n_local; i++) {
            for (int k = 0; k < 3; k++)
                if (local[i].pos[k] < mymin[k] || local[i].pos[k] > mymax[k]) inside = 0;
            if (domain_owner_pos(&d, local[i].pos) != rank) inside = 0;
        }

        int flags[3] = {disjoint, covers, inside}, all[3];
        MPI_Allreduce(flags, all, 3, MPI_INT, MPI_MIN, comm);
        ok_orb = all[0] && all[1] && all[2];

        if (rank == 0) {
            printf("  ORB    cajas disjuntas: %s\n", all[0] ? "si" : "NO");
            printf("  ORB    cubren el dominio: %s (vol %.6e vs %.6e)\n",
                   all[1] ? "si" : "NO", vsum, vglobal);
            printf("  ORB    particulas dentro de su caja: %s\n", all[2] ? "si" : "NO");
        }

        free(local);
        domain_free(&d);
        mpi_particle_type_free(&ptype);
    }

    /* --- b) Morton: el camino legado sigue siendo valido (se usa de linea base
             en los experimentos comparativos del informe) --- */
    {
        Particle *local = NULL;
        Domain d;
        MPI_Datatype ptype;
        mpi_particle_type_init(&ptype);
        domain_init(&d, comm);
        int n_local = setup_distributed_ex(&local, &capacity, &d, ptype, N, 42, 0, comm);

        int monotone = 1;
        for (int k = 0; k < P; k++)
            if (d.splitters[k] > d.splitters[k+1]) monotone = 0;
        int covers = (d.splitters[0] == 0 && d.splitters[P] == UINT64_MAX);
        int inside = 1;
        for (int i = 0; i < n_local; i++)
            if (local[i].morton < d.splitters[rank] ||
                local[i].morton >= d.splitters[rank+1]) inside = 0;

        int flags[3] = {monotone, covers, inside}, all[3];
        MPI_Allreduce(flags, all, 3, MPI_INT, MPI_MIN, comm);
        ok_morton = all[0] && all[1] && all[2];

        if (rank == 0)
            printf("  Morton splitters monotonos, cubrientes y respetados: %s\n",
                   ok_morton ? "si" : "NO");

        free(local);
        domain_free(&d);
        mpi_particle_type_free(&ptype);
    }

    int pass = ok_orb && ok_morton;
    if (rank == 0) printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
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

/* ------------------------------------------------------------------ */
/* Semana 4: LET y balance dinamico                                    */
/* ------------------------------------------------------------------ */

/*
 * Calculo de fuerzas por LET: arbol local sobre la caja global compartida,
 * intercambio de nodos esenciales, fantasmas insertados en ese mismo arbol y
 * fuerzas solo para el rango local [0, n_local). Devuelve n_ghost.
 */
static int let_forces(Particle **local, int n_local, int *capacity,
                      const Domain *d, MPI_Datatype ptype,
                      double theta, double softening, MPI_Comm comm) {
    Octree *t = octree_build_box(*local, n_local, d->gmin, d->gmax);
    octree_compute_mass(t, *local);

    int n_ghost = let_exchange(local, n_local, capacity, t, d,
                               theta, softening, ptype, comm);

    octree_insert_particles(t, *local, n_local, n_local + n_ghost);
    octree_compute_mass(t, *local);
    compute_forces_bh_range(t, *local, n_local + n_ghost, 0, n_local,
                            theta, softening);
    octree_free(t);
    return n_ghost;
}

/* Un paso completo del bucle hibrido usando LET en lugar de replicacion. */
static int distributed_step_let(Particle **local, int n_local, int *capacity,
                                Domain *d, MPI_Datatype ptype, int N,
                                double dt, double theta, double softening,
                                int use_work, int repartition, MPI_Comm comm) {
    leapfrog_kick(*local, n_local, dt / 2.0);
    leapfrog_drift(*local, n_local, dt);

    if (repartition) {
        domain_global_bounds(d, *local, n_local, comm);
        double gmin[3], gmax[3];
        for (int k = 0; k < 3; k++) { gmin[k] = d->gmin[k]; gmax[k] = d->gmax[k]; }
        morton_sort(*local, n_local, gmin, gmax);
        if (d->use_orb) domain_partition_orb(d, *local, n_local, use_work, comm);
        else            domain_partition_ex(d, *local, n_local, N, use_work, comm);
    }
    n_local = migrate_particles(local, n_local, capacity, d, ptype, comm);

    let_forces(local, n_local, capacity, d, ptype, theta, softening, comm);
    leapfrog_kick(*local, n_local, dt / 2.0);
    return n_local;
}

/*
 * Error relativo maximo de acc entre el arreglo global reunido y una referencia,
 * emparejando por id porque la migracion permuta las particulas.
 */
static double accel_error_by_id(const Particle *all, const Particle *ref,
                                const int *pos_of_id, int n) {
    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        const Particle *r = &ref[pos_of_id[all[i].id]];
        double diff[3];
        VEC3_SUB(diff, all[i].acc, r->acc);
        double den = VEC3_NORM(r->acc);
        if (den > 0.0) {
            double err = VEC3_NORM(diff) / den;
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

int test_let_exact_theta0(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 14: LET con theta=0 vs Barnes-Hut secuencial (bit a bit) ===\n");

    int N = 2000, capacity;
    double theta = 0.0, softening = 0.01;

    /* Referencia secuencial con el mismo theta=0: sin aproximaciones, el arbol
       se recorre entero hasta las hojas. */
    Particle *ref = malloc(N * sizeof(Particle));
    init_plummer(ref, N, 42);
    double mn[3], mx[3];
    octree_bounds(ref, N, mn, mx);
    morton_sort(ref, N, mn, mx);
    Octree *rt = octree_build(ref, N);
    octree_compute_mass(rt, ref);
    compute_forces_bh(rt, ref, N, theta, softening);
    octree_free(rt);

    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);
    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);

    int n_ghost = let_forces(&local, n_local, &capacity, &d, ptype,
                             theta, softening, comm);

    /* Con theta=0 no se acepta ningun resumen: cada proceso debe terminar
       viendo las N particulas del sistema. */
    int64_t nl = n_local, seen = n_local + n_ghost, min_seen;
    MPI_Allreduce(&seen, &min_seen, 1, MPI_INT64_T, MPI_MIN, comm);
    (void)nl;

    Particle *all = malloc(N * sizeof(Particle));
    int offset;
    replicate_particles(local, n_local, all, N, &offset, ptype, comm);

    int *pos_of_id = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) pos_of_id[ref[i].id] = i;
    double max_err = accel_error_by_id(all, ref, pos_of_id, N);

    int pass = (max_err == 0.0) && (min_seen == N);
    if (rank == 0) {
        printf("  Particulas vistas por proceso (locales+fantasmas): %lld de %d\n",
               (long long)min_seen, N);
        printf("  Error relativo maximo: %.6e\n", max_err);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(pos_of_id); free(all); free(local); free(ref);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_let_accuracy(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 15: Precision del LET a theta=0.5 ===\n");

    int N = 2000, capacity;
    double theta = 0.5, softening = 0.01;

    /* Dos referencias sobre el mismo estado: el metodo directo O(N^2), que es
       la verdad, y el Barnes-Hut secuencial, que es el error que ya se acepta. */
    Particle *exact = malloc(N * sizeof(Particle));
    init_plummer(exact, N, 42);
    double mn[3], mx[3];
    octree_bounds(exact, N, mn, mx);
    morton_sort(exact, N, mn, mx);

    Particle *seq = malloc(N * sizeof(Particle));
    memcpy(seq, exact, N * sizeof(Particle));

    compute_forces_direct(exact, N, softening);
    Octree *st = octree_build(seq, N);
    octree_compute_mass(st, seq);
    compute_forces_bh(st, seq, N, theta, softening);
    octree_free(st);

    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);
    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);
    let_forces(&local, n_local, &capacity, &d, ptype, theta, softening, comm);

    Particle *all = malloc(N * sizeof(Particle));
    int offset;
    replicate_particles(local, n_local, all, N, &offset, ptype, comm);

    int *pos_of_id = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) pos_of_id[exact[i].id] = i;

    double err_let_vs_seq = accel_error_by_id(all, seq,   pos_of_id, N);
    double err_let        = accel_error_by_id(all, exact, pos_of_id, N);

    /* El Barnes-Hut secuencial contra el directo: la vara con la que se mide. */
    double err_seq = 0.0;
    for (int i = 0; i < N; i++) {
        double diff[3];
        VEC3_SUB(diff, seq[i].acc, exact[i].acc);
        double den = VEC3_NORM(exact[i].acc);
        if (den > 0.0) {
            double e = VEC3_NORM(diff) / den;
            if (e > err_seq) err_seq = e;
        }
    }

    int pass = (err_let_vs_seq < 1e-3) && (err_let < 1.5 * err_seq);
    if (rank == 0) {
        printf("  LET vs Barnes-Hut secuencial: %.6e  (umbral 1e-3)\n", err_let_vs_seq);
        printf("  Barnes-Hut secuencial vs directo O(N^2): %.6e\n", err_seq);
        printf("  LET vs directo O(N^2):                   %.6e  (umbral %.6e)\n",
               err_let, 1.5 * err_seq);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(pos_of_id); free(all); free(local); free(seq); free(exact);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_let_volume(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);
    if (rank == 0)
        printf("=== Test 16: Volumen del LET (crecimiento de n_ghost con N) ===\n");

    const int Ns[2] = {5000, 20000};   /* factor 4 entre ambos */
    int64_t ghosts[2][2] = {{0,0},{0,0}};   /* [decomp][caso], 0=morton 1=orb */

    for (int dec = 0; dec < 2; dec++) {
        for (int c = 0; c < 2; c++) {
            int N = Ns[c], capacity;
            Particle *local = NULL;
            Domain d;
            MPI_Datatype ptype;
            mpi_particle_type_init(&ptype);
            domain_init(&d, comm);

            int n_local = setup_distributed_ex(&local, &capacity, &d, ptype,
                                               N, 42, dec, comm);
            int n_ghost = let_forces(&local, n_local, &capacity, &d, ptype,
                                     0.5, 0.01, comm);

            int64_t g = n_ghost, g_sum;
            MPI_Allreduce(&g, &g_sum, 1, MPI_INT64_T, MPI_SUM, comm);
            ghosts[dec][c] = g_sum / P;

            free(local);
            domain_free(&d);
            mpi_particle_type_free(&ptype);
        }
    }

    /* Si n_ghost creciera sublinealmente, el ratio fantasmas/locales bajaria al
       crecer N. Exponente: log(g2/g1) / log(N2/N1). 1.0 = lineal = no comprime. */
    double k_morton = log((double)ghosts[0][1] / (double)ghosts[0][0]) / log(4.0);
    double k_orb    = log((double)ghosts[1][1] / (double)ghosts[1][0]) / log(4.0);
    int pass = (k_orb < 0.9);

    if (rank == 0) {
        printf("  %-8s %10s %10s %10s %10s\n", "decomp", "g(N=5000)", "g(N=20000)",
               "ratio 5K", "exponente");
        printf("  %-8s %10lld %10lld %10.2f %10.2f\n", "morton",
               (long long)ghosts[0][0], (long long)ghosts[0][1],
               (double)ghosts[0][0] / (5000.0 / P), k_morton);
        printf("  %-8s %10lld %10lld %10.2f %10.2f\n", "orb",
               (long long)ghosts[1][0], (long long)ghosts[1][1],
               (double)ghosts[1][0] / (5000.0 / P), k_orb);
        printf("  Criterio: exponente ORB < 0.9 (semana 4 con Morton: 0.92)\n");
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }
    return pass;
}

int test_balance_work(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);
    if (rank == 0)
        printf("=== Test 17: Rebalanceo por costo vs por conteo ===\n");

    const int N = 20000, steps = 30;
    double imb[2] = {0, 0}, imb_n[2] = {0, 0};

    for (int mode = 0; mode < 2; mode++) {   /* 0 = conteo, 1 = costo */
        int capacity;
        Particle *local = NULL;
        Domain d;
        MPI_Datatype ptype;
        mpi_particle_type_init(&ptype);
        domain_init(&d, comm);

        int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);
        for (int s = 0; s < steps; s++)
            n_local = distributed_step_let(&local, n_local, &capacity, &d, ptype, N,
                                           0.001, 0.5, 0.01, mode,
                                           (s % 5 == 0), comm);

        /* Trabajo real del ultimo paso: determinista, a diferencia de los
           tiempos, que dependen de la maquina y del ruido del sistema. */
        int64_t w = 0;
        for (int i = 0; i < n_local; i++) w += local[i].work;

        int64_t w_max, w_sum, n_max, n_sum, nl = n_local;
        MPI_Allreduce(&w, &w_max, 1, MPI_INT64_T, MPI_MAX, comm);
        MPI_Allreduce(&w, &w_sum, 1, MPI_INT64_T, MPI_SUM, comm);
        MPI_Allreduce(&nl, &n_max, 1, MPI_INT64_T, MPI_MAX, comm);
        MPI_Allreduce(&nl, &n_sum, 1, MPI_INT64_T, MPI_SUM, comm);

        imb[mode]   = (double)w_max / ((double)w_sum / P);
        imb_n[mode] = (double)n_max / ((double)n_sum / P);

        free(local);
        domain_free(&d);
        mpi_particle_type_free(&ptype);
    }

    int pass = (imb[1] < imb[0]) && (imb[1] < 1.05);
    if (rank == 0) {
        printf("  Reparto por conteo: desbalance trabajo=%.4f  particulas=%.4f\n",
               imb[0], imb_n[0]);
        printf("  Reparto por costo:  desbalance trabajo=%.4f  particulas=%.4f\n",
               imb[1], imb_n[1]);
        printf("  (que el desbalance de PARTICULAS empeore es la senal esperada:\n");
        printf("   los procesos con zonas densas reciben menos particulas)\n");
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }
    return pass;
}

int test_let_energy_conservation(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
        printf("=== Test 18: Conservacion de energia con LET (200 pasos) ===\n");

    const int N = 5000, steps = 200;
    const double dt = 0.001, theta = 0.5, soft = 0.01;
    int capacity;
    Particle *local = NULL;
    Domain d;
    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);
    domain_init(&d, comm);

    int n_local = setup_distributed(&local, &capacity, &d, ptype, N, 42, comm);

    /* Energia inicial sobre el arbol del LET. */
    int n_ghost = let_forces(&local, n_local, &capacity, &d, ptype, theta, soft, comm);
    Octree *t = octree_build_box(local, n_local + n_ghost, d.gmin, d.gmax);
    octree_compute_mass(t, local);
    double e[2] = { kinetic_energy_range(local, n_local, 0, n_local),
                    potential_energy_bh_range(t, local, n_local + n_ghost,
                                              0, n_local, theta, soft) };
    octree_free(t);
    double g[2];
    MPI_Allreduce(e, g, 2, MPI_DOUBLE, MPI_SUM, comm);
    double E0 = g[0] + g[1];

    for (int s = 0; s < steps; s++)
        n_local = distributed_step_let(&local, n_local, &capacity, &d, ptype, N,
                                       dt, theta, soft, 1, (s % 20 == 0), comm);

    n_ghost = let_forces(&local, n_local, &capacity, &d, ptype, theta, soft, comm);
    t = octree_build_box(local, n_local + n_ghost, d.gmin, d.gmax);
    octree_compute_mass(t, local);
    e[0] = kinetic_energy_range(local, n_local, 0, n_local);
    e[1] = potential_energy_bh_range(t, local, n_local + n_ghost,
                                     0, n_local, theta, soft);
    octree_free(t);
    MPI_Allreduce(e, g, 2, MPI_DOUBLE, MPI_SUM, comm);
    double Ef = g[0] + g[1];

    double drift = fabs(Ef - E0) / fabs(E0);
    int pass = (drift < 0.01);
    if (rank == 0) {
        printf("  E0 = %.10e   Ef = %.10e\n", E0, Ef);
        printf("  Energy drift: %.6e  (umbral 1e-2)\n", drift);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }

    free(local);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return pass;
}

int test_let_volume_vs_procs(void) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);
    if (rank == 0)
        printf("=== Test 19: Volumen del LET, ORB vs Morton a igual configuracion ===\n");

    /* n_local constante (5000 por proceso) en vez de N constante. Es lo que
       aisla el efecto de la descomposicion: con N fijo, subir P achica n_local y
       el cociente fantasmas/locales crece por el denominador, no por la calidad
       del reparto. Medido con N fijo=20000, los fantasmas absolutos de ORB son
       casi constantes (4707/5901/5635 para P=2/4/8) y solo el cociente sube. */
    const int N = 5000 * P;
    int64_t ghosts[2];

    for (int dec = 0; dec < 2; dec++) {
        int capacity;
        Particle *local = NULL;
        Domain d;
        MPI_Datatype ptype;
        mpi_particle_type_init(&ptype);
        domain_init(&d, comm);

        int n_local = setup_distributed_ex(&local, &capacity, &d, ptype, N, 42, dec, comm);
        int n_ghost = let_forces(&local, n_local, &capacity, &d, ptype, 0.5, 0.01, comm);

        int64_t g = n_ghost, g_sum;
        MPI_Allreduce(&g, &g_sum, 1, MPI_INT64_T, MPI_SUM, comm);
        ghosts[dec] = g_sum / P;

        free(local);
        domain_free(&d);
        mpi_particle_type_free(&ptype);
    }

    /* La afirmacion que se verifica es comparativa y no absoluta: a igual
       configuracion, dominios compactos importan sustancialmente menos que
       tramos de curva. Medido entre 1.8x y 2.6x menos segun P. */
    double frac = (double)ghosts[1] / (double)ghosts[0];
    int pass = (frac < 0.6);

    if (rank == 0) {
        printf("  N=%d, P=%d, n_local=5000 por proceso\n", N, P);
        printf("  fantasmas/proceso con morton: %lld\n", (long long)ghosts[0]);
        printf("  fantasmas/proceso con orb:    %lld\n", (long long)ghosts[1]);
        printf("  ORB / Morton = %.2f  (criterio < 0.60)\n", frac);
        printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    }
    return pass;
}

#endif /* USE_MPI */

int run_all_tests(void) {
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int passed = 0, total = 19;

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
    collective += test_let_exact_theta0();
    collective += test_let_accuracy();
    collective += test_let_volume();
    collective += test_balance_work();
    collective += test_let_energy_conservation();
    collective += test_let_volume_vs_procs();

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
    printf("(los tests 10-19 requieren MPI: mpirun -np 4 ./nbody_mpi --validate)\n\n");
    return (passed == total) ? 0 : 1;
#endif
}
