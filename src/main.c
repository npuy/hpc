#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "force.h"
#include "leapfrog.h"
#include "morton.h"
#include "octree.h"
#include "barnes_hut.h"
#include "validation.h"
#include "vec3.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#include "mpi_types.h"
#include "domain.h"
#include "migration.h"
#include "metrics.h"
#include "let.h"
#endif

/* Parámetros de una corrida, comunes a la ruta secuencial y a la híbrida. */
typedef struct {
    int      N;
    double   dt;
    double   t_end;
    double   softening;
    char    *output;
    char    *init_type;
    int      use_bh;
    double   theta;
    unsigned seed;
    int      snapshot_interval;
    int      rebalance_every;
    int      show_metrics;
    int      use_let;       /* 1 = LET (etapa 8), 0 = replica global (semana 3) */
    int      balance_work;  /* 1 = reparto por costo, 0 = por conteo */
    double   theta_let;     /* theta del criterio de exportacion; <0 = usar theta */
} Config;

static void print_usage(const char *prog) {
    printf("Uso: %s [opciones]\n", prog);
    printf("  -n N          Numero de particulas (default: 1000)\n");
    printf("  -dt DT        Paso temporal (default: 0.001)\n");
    printf("  -t T_END      Tiempo final (default: 1.0)\n");
    printf("  -s EPS        Softening (default: 0.01)\n");
    printf("  -o FILE       Archivo de salida CSV\n");
    printf("  --init TYPE   Tipo de init: plummer|uniform|twobody (default: plummer)\n");
    printf("  --method M    Metodo de fuerza: direct|bh (default: bh)\n");
    printf("  --theta T     Angulo de apertura Barnes-Hut (default: 0.5)\n");
    printf("  --threads T   Hilos OpenMP (default: OMP_NUM_THREADS)\n");
    printf("  --rebalance-every K  Reparticion cada K pasos, solo MPI (default: 20)\n");
    printf("  --exchange E  Intercambio: let|replicate, solo MPI (default: let)\n");
    printf("  --balance B   Reparto: work|count, solo MPI (default: work)\n");
    printf("  --theta-let T Theta del criterio de exportacion LET (default: =theta)\n");
    printf("  --metrics     Desglose de tiempo por fase, solo MPI\n");
    printf("  --validate    Correr suite de validacion\n");
    printf("  --seed SEED   Semilla aleatoria (default: 42)\n");
}

/* Genera el estado inicial segun c->init_type. Puede ajustar *n (twobody). */
static void init_particles(Particle *p, int *n, const Config *c) {
    if (strcmp(c->init_type, "twobody") == 0)
        init_two_body(p, n);
    else if (strcmp(c->init_type, "uniform") == 0)
        init_uniform_cube(p, *n, c->seed);
    else
        init_plummer(p, *n, c->seed);
}

/* ------------------------------------------------------------------ */
/* Ruta secuencial (+ OpenMP)                                          */
/* ------------------------------------------------------------------ */

#ifndef USE_MPI

static int run_sequential(const Config *c) {
    int N = c->N;

    printf("N-Body Simulacion Secuencial\n");
    printf("  N = %d, dt = %g, t_end = %g, softening = %g\n",
           N, c->dt, c->t_end, c->softening);
    printf("  Init: %s, seed: %u\n", c->init_type, c->seed);
    if (c->use_bh) printf("  Metodo: Barnes-Hut (theta = %g)\n", c->theta);
    else           printf("  Metodo: directo O(N^2)\n");
#ifdef _OPENMP
    printf("  Hilos OpenMP: %d\n", omp_get_max_threads());
#endif
    printf("\n");

    Particle *particles = malloc(N * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Error: no se pudo asignar memoria para %d particulas\n", N);
        return 1;
    }
    init_particles(particles, &N, c);

    /* Morton sort inicial */
    double mn[3], mx[3];
    octree_bounds(particles, N, mn, mx);
    morton_sort(particles, N, mn, mx);

    double E0;
    if (c->use_bh) {
        Octree *tree = octree_build(particles, N);
        octree_compute_mass(tree, particles);
        compute_forces_bh(tree, particles, N, c->theta, c->softening);
        E0 = kinetic_energy(particles, N)
           + potential_energy_bh(tree, particles, N, c->theta, c->softening);
        octree_free(tree);
    } else {
        compute_forces_direct(particles, N, c->softening);
        E0 = kinetic_energy(particles, N) + potential_energy(particles, N, c->softening);
    }
    printf("Energia inicial E0 = %.10e\n", E0);

    if (c->output)
        write_particles_csv(c->output, particles, N, 0.0);

    struct timespec ts_start, ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int step = 0;
    int total_steps = (int)(c->t_end / c->dt);
    double Ef = E0;

    for (double t = 0; t < c->t_end - c->dt * 0.5; t += c->dt, step++) {
        leapfrog_kick(particles, N, c->dt / 2.0);
        leapfrog_drift(particles, N, c->dt);

        double E = 0.0;
        int is_snapshot = (c->snapshot_interval > 0 && step % c->snapshot_interval == 0);

        if (c->use_bh) {
            Octree *tree = octree_build(particles, N);
            octree_compute_mass(tree, particles);
            compute_forces_bh(tree, particles, N, c->theta, c->softening);
            if (is_snapshot)
                E = kinetic_energy(particles, N)
                  + potential_energy_bh(tree, particles, N, c->theta, c->softening);
            octree_free(tree);
        } else {
            compute_forces_direct(particles, N, c->softening);
            if (is_snapshot)
                E = kinetic_energy(particles, N)
                  + potential_energy(particles, N, c->softening);
        }

        leapfrog_kick(particles, N, c->dt / 2.0);

        if (is_snapshot) {
            double drift = fabs(E - E0) / fabs(E0);
            printf("  paso %d/%d  t=%.4f  E=%.10e  drift=%.4e\n",
                   step, total_steps, t + c->dt, E, drift);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double wall_time = (ts_end.tv_sec - ts_start.tv_sec)
                     + (ts_end.tv_nsec - ts_start.tv_nsec) * 1e-9;

    if (c->use_bh) {
        Octree *tree = octree_build(particles, N);
        octree_compute_mass(tree, particles);
        Ef = kinetic_energy(particles, N)
           + potential_energy_bh(tree, particles, N, c->theta, c->softening);
        octree_free(tree);
    } else {
        Ef = kinetic_energy(particles, N) + potential_energy(particles, N, c->softening);
    }

    printf("\n=== Resultado final ===\n");
    printf("  Pasos: %d\n", step);
    printf("  Energia final Ef = %.10e\n", Ef);
    printf("  Energy drift: %.6e\n", fabs(Ef - E0) / fabs(E0));
    printf("  Tiempo de ejecucion: %.3f s\n", wall_time);
    printf("  Tiempo por paso: %.3f ms\n", (wall_time / step) * 1000.0);

    if (c->output) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s_final.csv", c->output);
        write_particles_csv(fname, particles, N, c->t_end);
        printf("  Snapshot final: %s\n", fname);
    }

    free(particles);
    return 0;
}

#endif /* !USE_MPI */

/* ------------------------------------------------------------------ */
/* Ruta híbrida MPI + OpenMP                                           */
/* ------------------------------------------------------------------ */

#ifdef USE_MPI

/*
 * Recalcula el bounding box global, reordena localmente y recalcula los
 * splitters. Bounding box y splitters se actualizan SIEMPRE juntos: las claves
 * solo son comparables entre procesos si todos usan el mismo box, y los
 * splitters solo tienen sentido sobre el espacio de claves que las produjo.
 */
static void repartition(Domain *d, Particle *local, int n_local,
                        int64_t n_global, int use_work, MPI_Comm comm) {
    domain_global_bounds(d, local, n_local, comm);
    double gmin[3], gmax[3];
    for (int k = 0; k < 3; k++) { gmin[k] = d->gmin[k]; gmax[k] = d->gmax[k]; }
    morton_sort(local, n_local, gmin, gmax);
    domain_partition_ex(d, local, n_local, n_global, use_work, comm);
}

/*
 * Arma el arbol sobre el que se calculan las fuerzas del paso y devuelve
 * cuantas particulas contiene. Encapsula la diferencia entre los dos esquemas
 * de intercambio para que el bucle temporal no se bifurque:
 *
 *  - LET: el arbol se construye sobre las n_local locales con la caja global
 *    compartida, se exportan/importan los nodos esenciales y los fantasmas se
 *    insertan en ese mismo arbol (sin reconstruirlo). Las locales quedan en
 *    [0, n_local) y las fuerzas se calculan sobre ese rango.
 *
 *  - Replicacion: el arbol se construye sobre las N globales y las locales
 *    quedan en [offset, offset + n_local). Es el camino de la semana 3, que se
 *    conserva porque es el oraculo de los tests 15-16 y la linea base de todos
 *    los graficos comparativos del informe.
 *
 * Escribe en *base el arreglo sobre el que quedo armado el arbol y en *i0 el
 * indice donde empiezan las particulas locales dentro de el.
 */
static Octree *build_step_tree(Particle **local, int n_local, int *capacity,
                               Domain *d, const Config *c, MPI_Datatype ptype,
                               Particle *all, int N, Metrics *m,
                               Particle **base, int *i0, int *n_tree,
                               MPI_Comm comm) {
    double tmark;
    double theta_let = (c->theta_let >= 0.0) ? c->theta_let : c->theta;

    if (!c->use_let) {
        tmark = metrics_now();
        int offset;
        replicate_particles(*local, n_local, all, N, &offset, ptype, comm);
        m->comm_time += metrics_now() - tmark;

        tmark = metrics_now();
        Octree *t = octree_build(all, N);
        octree_compute_mass(t, all);
        m->tree_time += metrics_now() - tmark;

        *base = all; *i0 = offset; *n_tree = N;
        return t;
    }

    /* El arbol local se necesita completo (con masas y centros de masa) ANTES
       de exportar: el criterio de apertura conservador los usa. */
    tmark = metrics_now();
    Octree *t = octree_build_box(*local, n_local, d->gmin, d->gmax);
    octree_compute_mass(t, *local);
    m->tree_time += metrics_now() - tmark;

    tmark = metrics_now();
    int n_ghost = let_exchange(local, n_local, capacity, t, d,
                               theta_let, c->softening, ptype, comm);
    m->let_time += metrics_now() - tmark;
    m->ghost_sum += n_ghost;

    /* Insertar los fantasmas en el arbol ya construido en vez de rearmar uno
       fusionado: la caja raiz es la misma, asi que la estructura sigue valida y
       solo se subdividen las hojas donde caen los nuevos puntos. Las masas si
       hay que recalcularlas. */
    tmark = metrics_now();
    octree_insert_particles(t, *local, n_local, n_local + n_ghost);
    octree_compute_mass(t, *local);
    m->tree_time += metrics_now() - tmark;

    *base = *local; *i0 = 0; *n_tree = n_local + n_ghost;
    return t;
}

static int run_mpi(const Config *c) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    MPI_Datatype ptype;
    mpi_particle_type_init(&ptype);

    /* Rango 0 genera el estado inicial; N puede cambiar (twobody). */
    int N = c->N;
    Particle *gen = NULL;
    if (rank == 0) {
        gen = malloc(N * sizeof(Particle));
        init_particles(gen, &N, c);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, comm);

    if (rank == 0) {
        printf("N-Body Simulacion Hibrida MPI + OpenMP\n");
        printf("  N = %d, dt = %g, t_end = %g, softening = %g\n",
               N, c->dt, c->t_end, c->softening);
        printf("  Init: %s, seed: %u\n", c->init_type, c->seed);
        printf("  Metodo: Barnes-Hut (theta = %g)\n", c->theta);
        printf("  Procesos MPI: %d\n", P);
#ifdef _OPENMP
        printf("  Hilos OpenMP por proceso: %d\n", omp_get_max_threads());
#endif
        printf("  Intercambio: %s\n", c->use_let ? "LET" : "replicacion global");
        if (c->use_let && c->theta_let >= 0.0)
            printf("  Theta de exportacion LET: %g\n", c->theta_let);
        printf("  Reparto: por %s\n", c->balance_work ? "costo (work)" : "conteo");
        printf("  Reparticion cada %d pasos\n\n", c->rebalance_every);
    }

    /* Reparto inicial en bloques iguales (todavia sin criterio espacial). */
    int *counts = malloc(P * sizeof(int));
    int *displs = malloc(P * sizeof(int));
    for (int k = 0; k < P; k++) counts[k] = N / P + (k < N % P ? 1 : 0);
    displs[0] = 0;
    for (int k = 1; k < P; k++) displs[k] = displs[k-1] + counts[k-1];

    int n_local  = counts[rank];
    int capacity = (n_local > 0 ? n_local : 1);
    Particle *local = malloc(capacity * sizeof(Particle));
    MPI_Scatterv(gen, counts, displs, ptype, local, n_local, ptype, 0, comm);
    free(gen);
    free(counts);
    free(displs);

    /* Descomposicion espacial: cada proceso queda con un tramo Morton contiguo. */
    Domain d;
    domain_init(&d, comm);
    repartition(&d, local, n_local, N, 0, comm);
    n_local = migrate_particles(&local, n_local, &capacity, &d, ptype, comm);

    int64_t id_sum0, id_sum20;
    particle_checksum(local, n_local, &id_sum0, &id_sum20, comm);

    /* Buffer de replicacion: O(N) por proceso, el costo de no tener LET. En
       modo LET no se reserva — el ahorro de memoria es parte del punto. */
    Particle *all = c->use_let ? NULL : malloc(N * sizeof(Particle));

    Metrics m;
    metrics_reset(&m);

    /* Fuerzas y energia iniciales. */
    Particle *base; int i0, n_tree;
    Octree *tree = build_step_tree(&local, n_local, &capacity, &d, c, ptype,
                                   all, N, &m, &base, &i0, &n_tree, comm);
    compute_forces_bh_range(tree, base, n_tree, i0, i0 + n_local,
                            c->theta, c->softening);
    if (!c->use_let)
        for (int i = 0; i < n_local; i++) VEC3_COPY(local[i].acc, base[i0+i].acc);

    double e_local[2] = {
        kinetic_energy_range(local, n_local, 0, n_local),
        potential_energy_bh_range(tree, base, n_tree, i0, i0 + n_local,
                                  c->theta, c->softening)
    };
    octree_free(tree);
    double e_global[2];
    MPI_Allreduce(e_local, e_global, 2, MPI_DOUBLE, MPI_SUM, comm);
    double E0 = e_global[0] + e_global[1];
    if (rank == 0) printf("Energia inicial E0 = %.10e\n", E0);

    metrics_reset(&m);
    MPI_Barrier(comm);
    double t_start = metrics_now();

    int step = 0;
    int total_steps = (int)(c->t_end / c->dt);
    double Ef = E0, tmark;

    for (double t = 0; t < c->t_end - c->dt * 0.5; t += c->dt, step++) {
        int is_snapshot = (c->snapshot_interval > 0 && step % c->snapshot_interval == 0);

        tmark = metrics_now();
        leapfrog_kick(local, n_local, c->dt / 2.0);
        leapfrog_drift(local, n_local, c->dt);
        m.integrate_time += metrics_now() - tmark;

        /* Las particulas se movieron: reubicarlas antes de construir el arbol.
           Migrar al final del paso romperia la contiguidad de [offset, offset+n). */
        tmark = metrics_now();
        if (c->rebalance_every > 0 && step % c->rebalance_every == 0)
            repartition(&d, local, n_local, N, c->balance_work, comm);
        n_local = migrate_particles(&local, n_local, &capacity, &d, ptype, comm);
        m.migrate_time += metrics_now() - tmark;

        tree = build_step_tree(&local, n_local, &capacity, &d, c, ptype,
                               all, N, &m, &base, &i0, &n_tree, comm);

        tmark = metrics_now();
        compute_forces_bh_range(tree, base, n_tree, i0, i0 + n_local,
                                c->theta, c->softening);
        m.force_time += metrics_now() - tmark;

        if (!c->use_let)
            for (int i = 0; i < n_local; i++) VEC3_COPY(local[i].acc, base[i0+i].acc);

        double E = 0.0;
        if (is_snapshot) {
            e_local[0] = kinetic_energy_range(local, n_local, 0, n_local);
            e_local[1] = potential_energy_bh_range(tree, base, n_tree, i0, i0 + n_local,
                                                   c->theta, c->softening);
            MPI_Allreduce(e_local, e_global, 2, MPI_DOUBLE, MPI_SUM, comm);
            E = e_global[0] + e_global[1];
        }
        octree_free(tree);

        tmark = metrics_now();
        leapfrog_kick(local, n_local, c->dt / 2.0);
        m.integrate_time += metrics_now() - tmark;

        if (is_snapshot && rank == 0) {
            printf("  paso %d/%d  t=%.4f  E=%.10e  drift=%.4e\n",
                   step, total_steps, t + c->dt, E, fabs(E - E0) / fabs(E0));
        }
    }

    MPI_Barrier(comm);
    m.total_time = metrics_now() - t_start;
    m.n_local = n_local;
    m.steps = step;

    /* Energia final */
    Metrics discard;
    metrics_reset(&discard);
    tree = build_step_tree(&local, n_local, &capacity, &d, c, ptype,
                           all, N, &discard, &base, &i0, &n_tree, comm);
    e_local[0] = kinetic_energy_range(local, n_local, 0, n_local);
    e_local[1] = potential_energy_bh_range(tree, base, n_tree, i0, i0 + n_local,
                                           c->theta, c->softening);
    double mass_err = c->use_let
        ? let_mass_error(local, n_local, n_tree - n_local, comm) : 0.0;
    octree_free(tree);
    MPI_Allreduce(e_local, e_global, 2, MPI_DOUBLE, MPI_SUM, comm);
    Ef = e_global[0] + e_global[1];

    /* Conservacion de particulas e identidad tras toda la simulacion. */
    int64_t n_total;
    int64_t n_loc64 = n_local;
    MPI_Allreduce(&n_loc64, &n_total, 1, MPI_INT64_T, MPI_SUM, comm);
    int64_t id_sum, id_sum2;
    particle_checksum(local, n_local, &id_sum, &id_sum2, comm);

    if (rank == 0) {
        printf("\n=== Resultado final ===\n");
        printf("  Pasos: %d\n", step);
        printf("  Energia final Ef = %.10e\n", Ef);
        printf("  Energy drift: %.6e\n", fabs(Ef - E0) / fabs(E0));
        printf("  Tiempo de ejecucion: %.3f s\n", m.total_time);
        printf("  Tiempo por paso: %.3f ms\n", (m.total_time / step) * 1000.0);
        printf("  Particulas: %lld / %d  (%s)\n", (long long)n_total, N,
               n_total == N ? "conservadas" : "PERDIDAS");
        printf("  Checksum id: %s\n",
               (id_sum == id_sum0 && id_sum2 == id_sum20) ? "OK" : "ALTERADO");
        /* El LET resume pero no descarta: cada proceso debe "ver" la masa total
           del sistema entre sus locales y sus fantasmas. Es la verificacion mas
           barata de que la seleccion no se comio nada. */
        if (c->use_let)
            printf("  Error de masa del LET: %.3e\n", mass_err);
    }

    if (c->show_metrics) metrics_report(&m, comm);

    if (c->output) {
        /* La salida CSV necesita el estado global reunido en un proceso. En modo
           LET no hay buffer replicado durante la simulacion, asi que se arma uno
           solo para esto. */
        Particle *out = all ? all : malloc((size_t)N * sizeof(Particle));
        int off;
        replicate_particles(local, n_local, out, N, &off, ptype, comm);
        if (rank == 0) {
            char fname[256];
            snprintf(fname, sizeof(fname), "%s_final.csv", c->output);
            write_particles_csv(fname, out, N, c->t_end);
            printf("  Snapshot final: %s\n", fname);
        }
        if (!all) free(out);
    }

    free(all);
    free(local);
    domain_free(&d);
    mpi_particle_type_free(&ptype);
    return 0;
}

#endif /* USE_MPI */

/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    Config c = {
        .N = 1000, .dt = 0.001, .t_end = 1.0, .softening = 0.01,
        .output = NULL, .init_type = "plummer", .use_bh = 1, .theta = 0.5,
        .seed = 42, .snapshot_interval = 100, .rebalance_every = 20,
        .show_metrics = 0, .use_let = 1, .balance_work = 1, .theta_let = -1.0
    };
    int validate = 0;
    int threads = 0;
    char *method = "bh";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i+1 < argc)
            c.N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-dt") == 0 && i+1 < argc)
            c.dt = atof(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i+1 < argc)
            c.t_end = atof(argv[++i]);
        else if (strcmp(argv[i], "-s") == 0 && i+1 < argc)
            c.softening = atof(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc)
            c.output = argv[++i];
        else if (strcmp(argv[i], "--init") == 0 && i+1 < argc)
            c.init_type = argv[++i];
        else if (strcmp(argv[i], "--method") == 0 && i+1 < argc)
            method = argv[++i];
        else if (strcmp(argv[i], "--theta") == 0 && i+1 < argc)
            c.theta = atof(argv[++i]);
        else if (strcmp(argv[i], "--threads") == 0 && i+1 < argc)
            threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "--rebalance-every") == 0 && i+1 < argc)
            c.rebalance_every = atoi(argv[++i]);
        else if (strcmp(argv[i], "--exchange") == 0 && i+1 < argc)
            c.use_let = (strcmp(argv[++i], "replicate") != 0);
        else if (strcmp(argv[i], "--balance") == 0 && i+1 < argc)
            c.balance_work = (strcmp(argv[++i], "count") != 0);
        else if (strcmp(argv[i], "--theta-let") == 0 && i+1 < argc)
            c.theta_let = atof(argv[++i]);
        else if (strcmp(argv[i], "--metrics") == 0)
            c.show_metrics = 1;
        else if (strcmp(argv[i], "--validate") == 0)
            validate = 1;
        else if (strcmp(argv[i], "--seed") == 0 && i+1 < argc)
            c.seed = (unsigned)atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    c.use_bh = (strcmp(method, "bh") == 0);

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
#else
    (void)threads;
#endif

#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int rc;
    if (validate) {
        rc = run_all_tests();
    } else {
        if (!c.use_bh && rank == 0) {
            fprintf(stderr, "Aviso: --method direct no esta soportado en la version "
                            "distribuida (su rol es ser oraculo secuencial). "
                            "Usando Barnes-Hut.\n");
        }
        c.use_bh = 1;
        rc = run_mpi(&c);
    }
    MPI_Finalize();
    return rc;
#else
    if (validate) return run_all_tests();
    return run_sequential(&c);
#endif
}
