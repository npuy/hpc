# HPC - Facultad de Ingenieria (fing)

## Base de conocimiento

- [Configuracion y ejecucion de MPI](docs/mpi-setup.md) — Setup de mpich, compilacion con mpicc, ejecucion con mpirun, fix FI_PROVIDER=tcp, prerequisitos SSH para distribuido
- [Plan general de implementacion](docs/plan.md) — Arquitectura hibrida MPI/OpenMP, Barnes-Hut, etapas, cronograma 4 semanas, MVP, riesgos
- [Semana 1 — Plan detallado](docs/week1-plan.md) — Version secuencial: particulas, leapfrog KDK, fuerza O(N^2), Morton 3D, validacion, estructura de archivos
- [Semana 1 — Reporte](docs/week1-report.md) — Resultados: 4/4 tests validacion PASS, rendimiento O(N^2), componentes implementados, criterios de aceptacion
- [Semana 2 — Plan detallado](docs/week2-plan.md) — Barnes-Hut secuencial: octree (array-pool), centros de masa, criterio de apertura theta/MAC, fuerza O(N log N), validacion vs O(N^2)
- [Semana 2 — Reporte](docs/week2-report.md) — Resultados: 8/8 tests PASS, speedup BH vs directo (14.8x a N=100K), convergencia theta, conservacion energia, sin fugas
