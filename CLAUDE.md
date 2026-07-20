# HPC - Facultad de Ingenieria (fing)

## Base de conocimiento

- [Configuracion y ejecucion de MPI](docs/mpi-setup.md) — Setup de mpich, compilacion con mpicc, ejecucion con mpirun, fix FI_PROVIDER=tcp, prerequisitos SSH para distribuido
- [Plan general de implementacion](docs/plan.md) — Arquitectura hibrida MPI/OpenMP, Barnes-Hut, etapas, cronograma 4 semanas, MVP, riesgos
- [Semana 1 — Plan detallado](docs/week1-plan.md) — Version secuencial: particulas, leapfrog KDK, fuerza O(N^2), Morton 3D, validacion, estructura de archivos
- [Semana 1 — Reporte](docs/week1-report.md) — Resultados: 4/4 tests validacion PASS, rendimiento O(N^2), componentes implementados, criterios de aceptacion
- [Semana 2 — Plan detallado](docs/week2-plan.md) — Barnes-Hut secuencial: octree (array-pool), centros de masa, criterio de apertura theta/MAC, fuerza O(N log N), validacion vs O(N^2)
- [Semana 2 — Reporte](docs/week2-report.md) — Resultados: 8/8 tests PASS, speedup BH vs directo (14.8x a N=100K), convergencia theta, conservacion energia, sin fugas
- [Semana 3 — Plan detallado](docs/week3-plan.md) — Version hibrida MPI+OpenMP: replicacion via Allgatherv (LET queda para sem. 4), particionamiento Morton por biseccion de histograma, migracion con Alltoallv, tipo MPI de Particle, metricas por fase, tests 9-13
- [Semana 3 — Reporte](docs/week3-report.md) — Resultados: 13/13 tests PASS, fuerzas distribuidas bit a bit iguales a las secuenciales, 2 bugs corregidos (clamp Morton, carrera en fuerza directa). Hallazgo: el cuello de botella NO es el Allgatherv (1.3%) sino la construccion del arbol replicada por proceso
- [Spec: domain](docs/specs/domain.md) — Particionamiento por rangos Morton, biseccion de histograma simultanea, bounding box congelado, casos borde
- [Spec: migration](docs/specs/migration.md) — Alltoallv, orden de fases dentro del paso, replicacion provisional, checksums Sigma-id y Sigma-id2
- [Spec: metrics](docs/specs/metrics.md) — Timers por fase con MPI_Wtime, desbalance de particulas vs de trabajo
