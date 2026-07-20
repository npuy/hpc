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

    const char *names[6] = {"arbol", "fuerzas", "integracion",
                            "replica(Allgatherv)", "migracion(Alltoallv)", "total"};
    double local[6] = {m->tree_time, m->force_time, m->integrate_time,
                       m->comm_time, m->migrate_time, m->total_time};
    double vmin[6], vmax[6], vsum[6];

    MPI_Reduce(local, vmin, 6, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(local, vmax, 6, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(local, vsum, 6, MPI_DOUBLE, MPI_SUM, 0, comm);

    int64_t nmin, nmax, nsum;
    MPI_Reduce(&m->n_local, &nmin, 1, MPI_INT64_T, MPI_MIN, 0, comm);
    MPI_Reduce(&m->n_local, &nmax, 1, MPI_INT64_T, MPI_MAX, 0, comm);
    MPI_Reduce(&m->n_local, &nsum, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    if (rank != 0) return;

    printf("\n=== Metricas por fase (%d procesos, %lld pasos) ===\n",
           P, (long long)m->steps);
    printf("  %-22s %10s %10s %10s %8s\n", "fase", "min(s)", "max(s)", "avg(s)", "%total");
    double avg_total = vsum[5] / P;
    for (int i = 0; i < 6; i++) {
        double avg = vsum[i] / P;
        printf("  %-22s %10.4f %10.4f %10.4f %7.1f%%\n",
               names[i], vmin[i], vmax[i], avg,
               avg_total > 0 ? 100.0 * avg / avg_total : 0.0);
    }

    double avg_n = (double)nsum / P;
    double avg_force = vsum[1] / P;
    printf("\n  Particulas por proceso: min=%lld max=%lld avg=%.1f\n",
           (long long)nmin, (long long)nmax, avg_n);
    /* Desbalance de particulas: cuanto se desvia el proceso mas cargado del
       reparto ideal. Es lo que domain_partition busca mantener cerca de 1. */
    printf("  Desbalance de particulas (max/avg): %.4f\n",
           avg_n > 0 ? nmax / avg_n : 0.0);
    /* Desbalance de trabajo: mas informativo que el de particulas, porque el
       costo del recorrido no es uniforme entre regiones densas y vacias. */
    printf("  Desbalance de trabajo   (max/avg): %.4f\n",
           avg_force > 0 ? vmax[1] / avg_force : 0.0);
    /* Fraccion de comunicacion: el numero que motiva el LET de la semana 4. */
    printf("  Fraccion de comunicacion: %.1f%%\n",
           avg_total > 0 ? 100.0 * (vsum[3] + vsum[4]) / P / avg_total : 0.0);
}

#endif /* USE_MPI */
