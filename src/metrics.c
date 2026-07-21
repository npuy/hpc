#include "metrics.h"

#ifdef USE_MPI

#include <stdio.h>
#include <string.h>

void metrics_reset(Metrics *m) {
    memset(m, 0, sizeof(Metrics));
}

double metrics_now(void) {
    return MPI_Wtime();
}

void metrics_report(const Metrics *m, MPI_Comm comm) {
    int rank, P;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &P);

    enum { NPHASE = 7, I_FORCE = 1, I_TOTAL = 6 };
    const char *names[NPHASE] = {"arbol", "fuerzas", "integracion",
                                 "replica(Allgatherv)", "let(seleccion+Alltoallv)",
                                 "migracion(Alltoallv)", "total"};
    double local[NPHASE] = {m->tree_time, m->force_time, m->integrate_time,
                            m->comm_time, m->let_time, m->migrate_time,
                            m->total_time};
    double vmin[NPHASE], vmax[NPHASE], vsum[NPHASE];

    MPI_Reduce(local, vmin, NPHASE, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(local, vmax, NPHASE, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(local, vsum, NPHASE, MPI_DOUBLE, MPI_SUM, 0, comm);

    int64_t nmin, nmax, nsum, gsum, gmax;
    MPI_Reduce(&m->n_local, &nmin, 1, MPI_INT64_T, MPI_MIN, 0, comm);
    MPI_Reduce(&m->n_local, &nmax, 1, MPI_INT64_T, MPI_MAX, 0, comm);
    MPI_Reduce(&m->n_local, &nsum, 1, MPI_INT64_T, MPI_SUM, 0, comm);
    MPI_Reduce(&m->ghost_sum, &gsum, 1, MPI_INT64_T, MPI_SUM, 0, comm);
    MPI_Reduce(&m->ghost_sum, &gmax, 1, MPI_INT64_T, MPI_MAX, 0, comm);

    if (rank != 0) return;

    printf("\n=== Metricas por fase (%d procesos, %lld pasos) ===\n",
           P, (long long)m->steps);
    printf("  %-26s %10s %10s %10s %8s\n", "fase", "min(s)", "max(s)", "avg(s)", "%total");
    double avg_total = vsum[I_TOTAL] / P;
    for (int i = 0; i < NPHASE; i++) {
        double avg = vsum[i] / P;
        printf("  %-26s %10.4f %10.4f %10.4f %7.1f%%\n",
               names[i], vmin[i], vmax[i], avg,
               avg_total > 0 ? 100.0 * avg / avg_total : 0.0);
    }

    double avg_n = (double)nsum / P;
    double avg_force = vsum[I_FORCE] / P;

    /* Volumen del LET: cuántos fantasmas importó en promedio cada proceso por
       paso. Es la métrica que dice si el LET escala — debe crecer mucho más
       despacio que N, no proporcionalmente como la réplica. */
    if (m->steps > 0 && gsum > 0) {
        double avg_ghost = (double)gsum / P / (double)m->steps;
        printf("\n  Fantasmas por proceso y paso: avg=%.0f  max=%.0f\n",
               avg_ghost, (double)gmax / (double)m->steps);
        printf("  Fantasmas / particulas locales: %.2f\n",
               avg_n > 0 ? avg_ghost / avg_n : 0.0);
    }
    printf("  Particulas por proceso: min=%lld max=%lld avg=%.1f\n",
           (long long)nmin, (long long)nmax, avg_n);
    /* Desbalance de particulas: cuanto se desvia el proceso mas cargado del
       reparto ideal. Es lo que domain_partition busca mantener cerca de 1. */
    printf("  Desbalance de particulas (max/avg): %.4f\n",
           avg_n > 0 ? nmax / avg_n : 0.0);
    /* Desbalance de trabajo: mas informativo que el de particulas, porque el
       costo del recorrido no es uniforme entre regiones densas y vacias. */
    printf("  Desbalance de trabajo   (max/avg): %.4f\n",
           avg_force > 0 ? vmax[1] / avg_force : 0.0);
    /* Fraccion de comunicacion. La semana 3 midio 1.3% con replicacion: es el
       numero que descarto a la red como cuello de botella y redirigio el
       argumento del LET hacia el costo del arbol replicado. */
    printf("  Fraccion de comunicacion: %.1f%%\n",
           avg_total > 0 ? 100.0 * (vsum[3] + vsum[4] + vsum[5]) / P / avg_total : 0.0);
}

#endif /* USE_MPI */
