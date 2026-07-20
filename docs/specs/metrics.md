# Módulo: metrics

## Descripción general

Instrumenta el paso de simulación separándolo en fases y agrega los tiempos de todos los procesos. Es lo que convierte una corrida en un dato utilizable para el informe: sin el desglose, un speedup pobre no dice **dónde** está el problema.

## Conceptos

### Por qué el desglose por fase

El número que importa no es el tiempo total sino su composición. En este proyecto el desglose es lo que permitió identificar el cuello de botella real de la semana 3: la fracción de comunicación resultó baja (~1% el `Allgatherv`) mientras que la construcción del árbol —replicada en cada proceso— crecía de 0.21 s a 1.06 s al pasar de 1 a 8 procesos. Sin `tree_time` medido por separado, ese diagnóstico habría sido invisible y se habría atribuido erróneamente a la comunicación.

### `MPI_Wtime`, no `clock_gettime`

En el camino MPI todos los timers usan `MPI_Wtime()`: es la fuente de tiempo que MPI garantiza consistente entre procesos, mientras que `clock_gettime` puede diferir entre nodos de un cluster.

### `MPI_Barrier` antes de medir

Las fases colectivas se miden después de un `MPI_Barrier`. Sin la barrera, el tiempo de una comunicación colectiva absorbe el desbalance de cómputo previo: un proceso que llegó tarde al `Allgatherv` haría parecer que la comunicación es lenta, cuando en realidad estuvo esperando a otro. El barrier separa "cuánto tardó comunicar" de "cuánto esperé a los demás".

En el bucle principal la barrera se usa al abrir y cerrar la medición total; dentro del paso las fases se acumulan sin barrera para no penalizar la corrida, con la salvedad de que `migrate_time` incluye la espera en el colectivo (visible como una diferencia grande entre `min` y `max`).

### Las dos métricas de desbalance

| Métrica | Qué mide | Para qué |
|---|---|---|
| `max/avg` de `n_local` | Desbalance de **partículas** | Verifica que `domain_partition` reparta bien |
| `max/avg` de `force_time` | Desbalance de **trabajo** | El que realmente importa |

Son distintas y la segunda es la informativa: el costo del recorrido Barnes-Hut no es uniforme entre partículas, así que un reparto perfecto en cantidad (medido: 1.0000) puede seguir dando desbalance de trabajo (medido: 1.11). Esa brecha es el argumento cuantitativo para el rebalanceo por costo de la semana 4, que debería repartir por trabajo estimado y no por conteo.

## API

| Función | Rol |
|---|---|
| `metrics_reset` | Pone los contadores en cero |
| `metrics_now` | Marca de tiempo (`MPI_Wtime`) |
| `metrics_report` | Agrega MIN/MAX/AVG sobre los procesos e imprime en el rango 0 |

## Estructura `Metrics`

| Campo | Fase |
|---|---|
| `tree_time` | `octree_build` + `octree_compute_mass` |
| `force_time` | Recorrido Barnes-Hut |
| `integrate_time` | `leapfrog_kick` + `leapfrog_drift` |
| `comm_time` | `replicate_particles` (`Allgatherv`) |
| `migrate_time` | `migrate_particles` (`Alltoallv` + reordenamiento) |
| `total_time` | Pared, del bucle completo |
| `n_local`, `steps` | Contexto para los cocientes de desbalance |

`metrics_report` es **colectiva** (usa `MPI_Reduce`): todos los procesos deben invocarla, aunque solo el rango 0 imprima.

## Uso

```bash
mpirun -np 4 ./nbody_mpi -n 100000 -t 0.01 --metrics
```
