#ifndef METRICS_H
#define METRICS_H

/*
 * metrics.h — Instrumentación por fase de la simulación híbrida.
 *
 * Separa el tiempo del paso en cómputo (árbol, fuerzas, integración) y
 * comunicación (réplica, migración). El desglose es lo que permite argumentar
 * dónde está el cuello de botella: si comm_time crece con P mientras
 * force_time baja, el límite es la comunicación y no el cómputo — exactamente
 * el diagnóstico que motiva el LET de la semana 4.
 *
 * metrics_report agrega los valores de todos los procesos con MIN/MAX/AVG. El
 * cociente max/avg de force_time es la métrica de desbalance de carga que
 * justifica el rebalanceo dinámico.
 */

#ifdef USE_MPI

#include <mpi.h>
#include <stdint.h>

typedef struct {
    double  tree_time;       /* octree_build + octree_compute_mass */
    double  force_time;      /* recorrido Barnes-Hut */
    double  integrate_time;  /* kick + drift */
    double  comm_time;       /* replicate_particles (Allgatherv) */
    double  migrate_time;    /* migrate_particles (Alltoallv + reordenamiento) */
    double  total_time;      /* pared, del bucle completo */
    int64_t n_local;         /* partículas locales al final (para el desbalance) */
    int64_t steps;
} Metrics;

/* Pone todos los contadores en cero. */
void metrics_reset(Metrics *m);

/*
 * Marca de tiempo para acumular fases. Usa MPI_Wtime, que es la fuente de
 * tiempo consistente entre procesos.
 */
double metrics_now(void);

/*
 * Imprime (solo en el rango 0) el desglose por fase con min/max/promedio sobre
 * los procesos, el desbalance de carga y el porcentaje de tiempo en
 * comunicación. Colectivo: todos los procesos deben llamarlo.
 */
void metrics_report(const Metrics *m, MPI_Comm comm);

#endif /* USE_MPI */

#endif
