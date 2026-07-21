# HPC - Facultad de Ingenieria (fing)

## Base de conocimiento

- [Arquitectura general](docs/arquitectura.md) — Como funciona todo: paso de simulacion, las tres ideas (MAC, ORB, LET), mapa de modulos, struct Particle, caminos alternativos por flag, los 19 tests, CLI, resultados por semana
- [Guia de ejecucion en pcunix](docs/pcunix-guia.md) — Estado real del cluster (pcunix40 caida, 2 nucleos fisicos por nodo, home NFS compartido), colocacion de procesos, puesta a punto, que ya se verifico (19/19 con GCC+MPICH), como seria la corrida completa y que mirar
- [Runbook: ejecutar en pcunix paso a paso](docs/pcunix-runbook.md) — Guia operativa: 10 pasos con comando exacto, salida esperada y que hacer si falla; calibracion antes de lanzar la matriz; corrida por bloques; tabla de problemas conocidos
- [Configuracion y ejecucion de MPI](docs/mpi-setup.md) — Setup de mpich, compilacion con mpicc, ejecucion con mpirun, fix FI_PROVIDER=tcp, prerequisitos SSH para distribuido
- [REPORTE FINAL](docs/final-report.md) — Sintesis del proyecto completo, organizada por tema y no por semana: problema, solucion construida (nucleo secuencial, reparto hibrido, ORB, migracion, LET, balance), justificacion de las primitivas MPI y OpenMP usadas, metodologia de validacion (escalera de oraculos, tests bit a bit), los 6 hallazgos, resultados consolidados en pcunix, 8 lecciones metodologicas
- [Plan general de implementacion](docs/plan.md) — Arquitectura hibrida MPI/OpenMP, Barnes-Hut, etapas, cronograma 4 semanas, MVP, riesgos
- [Semana 1 — Plan detallado](docs/week1-plan.md) — Version secuencial: particulas, leapfrog KDK, fuerza O(N^2), Morton 3D, validacion, estructura de archivos
- [Semana 1 — Reporte](docs/week1-report.md) — Resultados: 4/4 tests validacion PASS, rendimiento O(N^2), componentes implementados, criterios de aceptacion
- [Semana 2 — Plan detallado](docs/week2-plan.md) — Barnes-Hut secuencial: octree (array-pool), centros de masa, criterio de apertura theta/MAC, fuerza O(N log N), validacion vs O(N^2)
- [Semana 2 — Reporte](docs/week2-report.md) — Resultados: 8/8 tests PASS, speedup BH vs directo (14.8x a N=100K), convergencia theta, conservacion energia, sin fugas
- [Semana 3 — Plan detallado](docs/week3-plan.md) — Version hibrida MPI+OpenMP: replicacion via Allgatherv (LET queda para sem. 4), particionamiento Morton por biseccion de histograma, migracion con Alltoallv, tipo MPI de Particle, metricas por fase, tests 9-13
- [Semana 3 — Reporte](docs/week3-report.md) — Resultados: 13/13 tests PASS, fuerzas distribuidas bit a bit iguales a las secuenciales, 2 bugs corregidos (clamp Morton, carrera en fuerza directa). Hallazgo: el cuello de botella NO es el Allgatherv (1.3%) sino la construccion del arbol replicada por proceso
- [Semana 4 — Plan detallado](docs/week4-plan.md) — LET (criterio de exportacion conservador, fantasmas como pseudo-particulas), octree_build_box con caja compartida, balance por costo (campo work + biseccion ponderada), fin de la exactitud bit a bit y nueva escalera de validacion, tests 14-18, experimentos e informe final
- [Semana 4 — Reporte](docs/week4-report.md) — Resultados: 17/18 tests, LET correcto (theta=0 bit a bit) pero NO comprime (n_ghost ~ N^0.92, lineal). Causa medida: un tramo Morton contiguo no es una region compacta, el AABB del destino cubre 44.6% del dominio. Balance por costo si cumple (desbalance de trabajo 1.21 -> 1.00). Arreglo pendiente: ORB
- [Semana 5 — Plan detallado](docs/week5-plan.md) — ORB (biseccion recursiva ortogonal) para que los dominios sean cajas y el LET comprima: arbol ORB replicado, migracion por posicion en vez de por clave, Morton degradado a localidad de cache, harness de medicion GO/NO-GO el dia 1, test 16 pasa a ser criterio de aceptacion (k<0.9), tests 12 reescrito y 19 nuevo
- [Semana 5 — Reporte](docs/week5-report.md) — ORB implementado: los dos criterios que fallaban ahora se cumplen (n_ghost ~ N^0.66 vs 0.92; tree_time crece 1.6x vs 4.8x). 19/19 tests con P=2,3,4,8. let.c no necesito cambios. Brecha Plummer/uniform de 2.77x a 1.83x, que confirma el diagnostico de la semana 4
- [Resultados medidos en pcunix](docs/resultados-pcunix.txt) — Salida cruda de la corrida completa (8 bloques, 6m40s, 5 nodos i3-4170): 19/19 con P=2,3,4,8 distribuidos; n_ghost k=0.67 orb vs 0.93 morton; speedup 3.20x con P=4
- [Spec: orb](docs/specs/orb.md) — Arbol ORB replicado, eje por dimension mas larga, biseccion del trabajo, pertenencia por camino recorrido y no por geometria, Morton degradado a localidad de cache
- [Spec: let](docs/specs/let.md) — Caja raiz compartida, MAC conservador por distancia minima a la caja destino, fantasmas como pseudo-particulas, orden de fases, limitacion conocida
- [Spec: balance](docs/specs/balance.md) — Campo work medido en el traversal, biseccion ponderada por prefijo, caso del primer paso, inversion esperada de los desbalances
- [Spec: domain](docs/specs/domain.md) — Particionamiento por rangos Morton, biseccion de histograma simultanea, bounding box congelado, casos borde
- [Spec: migration](docs/specs/migration.md) — Alltoallv, orden de fases dentro del paso, replicacion provisional, checksums Sigma-id y Sigma-id2
- [Spec: metrics](docs/specs/metrics.md) — Timers por fase con MPI_Wtime, desbalance de particulas vs de trabajo
