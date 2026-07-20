# Semana 3 — Reporte: Versión Híbrida MPI + OpenMP

## Archivos implementados

```
src/
  mpi_types.h / .c    -- tipo derivado de MPI para Particle (struct + resized)
  domain.h / .c       -- particionamiento Morton por biseccion de histograma
  migration.h / .c    -- migracion (Alltoallv), replicacion (Allgatherv), checksums
  metrics.h / .c      -- timers por fase, agregacion MIN/MAX/AVG entre procesos
  barnes_hut.h / .c   -- OpenMP + compute_forces_bh_range (composicion del hibrido)
  force.h / .c        -- compute_forces_direct_par (sin simetria, paralelizable)
  leapfrog.c          -- OpenMP en kick y drift
  morton.h / .c       -- clamp en morton_encode (bug), morton_keys separado
  particle.h / .c     -- campo id, CSV con columna id
  main.c              -- run_sequential / run_mpi bajo #ifdef USE_MPI
  validation.h / .c   -- +5 tests (9-13)
Makefile              -- deteccion de OpenMP en macOS/Linux, target nbody_mpi
docs/specs/
  domain.md, migration.md, metrics.md
```

## Decisión de arquitectura: replicación en lugar de LET

El LET es etapa 8 (semana 4). Para tener una versión híbrida funcionando sin depender de él, cada proceso **replica** el arreglo global vía `MPI_Allgatherv`, construye el árbol completo, pero **solo calcula fuerzas para su tramo** de partículas (`compute_forces_bh_range`).

La consecuencia buscada es que las fuerzas resultan **numéricamente idénticas** a la versión secuencial: el test 13 mide error relativo máximo **exactamente 0.0** (bit a bit). La consecuencia no buscada —y el hallazgo principal de la semana— se documenta más abajo.

## Dos bugs corregidos, ambos detectados al revisar el código de las semanas 1-2

### 1. `morton_encode` no acotaba las coordenadas (`src/morton.c`)

La coordenada normalizada se casteaba a `uint64_t` sin acotarla. Con el bounding box congelado entre reparticiones, una partícula que se sale del box producía:

- si excedía el máximo: `expand_bits` la enmascara con `& 0x1fffff` y **envuelve** el valor — una partícula del extremo derecho recibía una clave del extremo izquierdo y migraba al proceso equivocado;
- si quedaba por debajo del mínimo: el double negativo convertido a `uint64_t` es **comportamiento indefinido** en C.

No se manifestaba en la semana 2 porque el box se recalculaba en cada paso y siempre contenía todo. Corregido acotando antes del cast.

### 2. `compute_forces_direct` no es paralelizable (`src/force.c`)

Aprovecha la tercera ley de Newton acumulando en `p[j].acc` dentro del bucle interno. Con `parallel for` sobre `i`, dos hilos escriben el mismo `p[j].acc`: carrera de datos silenciosa.

Se agregó `compute_forces_direct_par`, sin simetría (cada `i` escribe solo `p[i].acc`), al doble de flops pero sin carreras y bit a bit reproducible. La versión simétrica queda intacta como oráculo secuencial de los tests 1-8.

## Resultados de validación

Suite completa: **13/13 PASS** (verificado con `-np 4` y `-np 8`).

| Test | Descripción | Criterio | Resultado | Valor obtenido |
|------|-------------|----------|-----------|----------------|
| 1-8 | Suite de semanas 1-2 | — | **PASS** | sin regresión |
| 9 | Determinismo OpenMP (1/2/4/8 hilos) | bit a bit | **PASS** | idéntico en `bh` y `directo-par` |
| 10 | Conservación de partículas (100 pasos) | Σ n_local == N | **PASS** | siempre 5000 |
| 11 | Identidad (sin pérdidas ni duplicados) | Σid, Σid² invariantes | **PASS** | 12497500 / 41654167500 |
| 12 | Particiones disjuntas y cubrientes | válidas + desbalance < 1.15 | **PASS** | desbalance **1.0000** |
| 13 | Fuerzas distribuidas vs secuenciales | error < 1e-12 | **PASS** | **0.0 exacto** |

El error exactamente nulo del test 13 confirma que el esquema de replicación reproduce el árbol secuencial idénticamente, y que el emparejamiento por `id` tras la permutación de la migración es correcto.

## Entorno de medición

> **Importante:** los benchmarks se corrieron en un portátil Apple Silicon con **4 núcleos performance + 4 efficiency** (8 lógicos), Open MPI 5.0.9, clang + libomp. **No** en las máquinas `pcunix` de fing. La heterogeneidad de núcleos afecta directamente las mediciones con más de 4 hilos/procesos y está señalada donde corresponde. Los números deben re-medirse en pcunix (núcleos homogéneos) antes de ir al informe final.

## a) Escalabilidad OpenMP (1 proceso, N=100K)

| Hilos | ms/paso | Speedup total | Fase fuerzas (s) | Speedup fuerzas | Fase árbol (s) |
|-------|---------|---------------|------------------|-----------------|----------------|
| 1 | 1410.0 | 1.00× | 12.604 | 1.00× | 0.230 |
| 2 | 737.0 | 1.91× | 6.586 | 1.91× | 0.243 |
| 4 | 409.7 | 3.44× | 3.594 | 3.51× | 0.236 |
| 8 | 358.1 | **3.94×** | 2.894 | **4.36×** | 0.236 |

**El criterio de aceptación (>5× con 8 hilos) no se cumple: 3.94×.** El desglose identifica la causa sin ambigüedad.

La fase de fuerzas es embarazosamente paralela y no tiene componente serial, pero aun así **solo alcanza 4.36× con 8 hilos**. Eso descarta a Amdahl como explicación y apunta al hardware: con 4 hilos (que caben en los 4 núcleos *performance*) el speedup es 3.51× — 88% de eficiencia, muy bueno —, mientras que agregar los 4 núcleos *efficiency* solo aporta 0.85× más, es decir cada núcleo E rinde ~21% de un núcleo P.

La parte serial existe y es medible: `arbol` se mantiene constante en ~0.236 s sin importar el número de hilos (7% del total a 8 hilos), como corresponde a `octree_build`, que no está paralelizado.

**Estimación para hardware homogéneo:** si la fase de fuerzas escalara linealmente, con 7% de fracción serial Amdahl daría ≈ 1/(0.07 + 0.93/8) ≈ **5.4×** a 8 hilos, que sí cumpliría el criterio. Es una proyección, no una medición: hay que verificarla en pcunix.

## b) Escalabilidad MPI (1 hilo, N=100K)

| Procesos | Total (s) | Speedup | Árbol (s) | Fuerzas (s) | Allgatherv (s) | Alltoallv (s) |
|----------|-----------|---------|-----------|-------------|----------------|---------------|
| 1 | 14.330 | 1.00× | 0.235 | 12.710 | 0.002 | 0.067 |
| 2 | 7.828 | 1.83× | 0.267 | 6.532 | 0.007 | 0.272 |
| 4 | 5.120 | **2.80×** | 0.667 | 3.655 | 0.017 | 0.315 |
| 8 | 5.480 | 2.61× | 1.133 | 3.270 | 0.071 | 0.538 |

**El criterio (>3× con 4 procesos) no se cumple: 2.80×.** Y con 8 procesos el rendimiento **empeora** respecto a 4.

La fase de fuerzas sí escala bien: 12.710 → 3.655 s con 4 procesos es **3.48×** (87% de eficiencia). Lo que arrastra el total hacia abajo es la columna `árbol`, que **crece 4.8× al pasar de 1 a 8 procesos** (0.235 → 1.133 s).

### El hallazgo principal: el cuello de botella no es la comunicación

El plan de la semana anticipaba que el `Allgatherv` dominaría y que ese sería el argumento para el LET. **La medición lo refuta.** El `Allgatherv` es barato: 0.071 s con 8 procesos, apenas **1.3%** del tiempo total.

El costo real de la replicación es la **construcción del árbol duplicada**. Cada proceso construye el árbol completo sobre las N partículas: trabajo idéntico, hecho P veces en paralelo, que no se reparte entre nadie y además satura el ancho de banda de memoria cuando P procesos lo hacen simultáneamente. Es trabajo redundante puro, y crece con P.

El LET resuelve ambas cosas, pero la razón de peso es esta segunda: con LET cada proceso construye un árbol sobre sus N/P partículas locales más los nodos importados, no sobre las N.

## c) Configuración híbrida (producto P×T = 8, N=100K)

| Config | Total (s) | Árbol (s) | Fuerzas (s) |
|--------|-----------|-----------|-------------|
| 1 proc × 8 hilos | **3.451** | 0.207 | 2.916 |
| 2 proc × 4 hilos | 3.457 | 0.233 | 2.815 |
| 4 proc × 2 hilos | 3.990 | 0.460 | 2.940 |
| 8 proc × 1 hilo | 5.192 | 1.064 | 3.156 |

La fase de fuerzas es prácticamente **constante** (2.8–3.2 s) en las cuatro configuraciones: el trabajo paralelo se reparte igual de bien por hilos que por procesos. Toda la diferencia está en `árbol`, que crece 5× de 1 a 8 procesos.

**Conclusión operativa: en el esquema replicado conviene minimizar procesos MPI y maximizar hilos OpenMP**, porque los hilos comparten un único árbol mientras que los procesos lo duplican. El óptimo medido es 1 proceso × 8 hilos (empatado con 2×4).

Esto es una propiedad del esquema provisional, no del algoritmo: es exactamente lo que el LET debe invertir en la semana 4.

## d) Escalabilidad débil (N/P = 25K constante, 1 hilo)

| Procesos | N | Total (s) | Eficiencia débil | Árbol (s) | Fuerzas (s) | Allgatherv (s) |
|----------|---|-----------|------------------|-----------|-------------|----------------|
| 1 | 25.000 | 2.819 | 100% | 0.048 | 2.503 | 0.000 |
| 2 | 50.000 | 3.532 | 80% | 0.187 | 2.901 | 0.003 |
| 4 | 100.000 | 4.856 | 58% | 0.286 | 3.757 | 0.016 |
| 8 | 200.000 | 11.780 | **24%** | 2.592 | 6.945 | 0.133 |

La degradación era esperada y su forma confirma el diagnóstico: con replicación, el trabajo por proceso **no** es constante al mantener N/P fijo. Cada proceso construye un árbol sobre las N totales, y N crece linealmente con P — de ahí que `árbol` pase de 0.048 s a 2.592 s (54×). La fase de fuerzas también crece, pero solo como `log N` (2.503 → 6.945 s), que es el comportamiento correcto de Barnes-Hut.

Es la demostración cuantitativa de que la replicación no es una estrategia escalable, solo un andamio de transición.

## Desglose por fase de una corrida representativa

N=10.000, 4 procesos × 2 hilos, 200 pasos:

| Fase | min (s) | max (s) | avg (s) | % total |
|------|---------|---------|---------|---------|
| árbol | 0.487 | 0.495 | 0.491 | 9.4% |
| fuerzas | 3.157 | 4.257 | 3.818 | 73.4% |
| integración | 0.016 | 0.034 | 0.024 | 0.5% |
| réplica (Allgatherv) | 0.031 | 0.037 | 0.033 | 0.6% |
| migración (Alltoallv) | 0.338 | 1.367 | 0.750 | 14.4% |

- Desbalance de **partículas** (max/avg): **1.0024**
- Desbalance de **trabajo** (max/avg): **1.1149**

Los dos desbalances difieren, y esa brecha es informativa: `domain_partition` reparte las partículas casi perfectamente, pero el costo del recorrido Barnes-Hut no es uniforme (una partícula en zona densa abre más nodos), así que el trabajo real queda 11% desbalanceado. **El rebalanceo de la semana 4 debería repartir por costo estimado, no por conteo de partículas.**

El amplio rango min/max de `migración` (0.34 a 1.37 s) es en buena parte espera dentro del colectivo: los procesos que terminan antes su cómputo se bloquean en el `Alltoallv` esperando a los demás.

## Conservación de energía y robustez

| Configuración | Pasos | Energy drift | Partículas | Checksum id |
|---|---|---|---|---|
| N=10.000, 4×2, dt=0.001 | 1000 | 6.51e-05 | 10000/10000 | OK |
| N=100.000, 4×1, dt=0.001 | 10 | 4.35e-06 | conservadas | OK |
| N=500.000, 8×1, dt=0.001 | 5 | 3.53e-06 | 500000/500000 | OK |

Casos borde verificados:

| Caso | Resultado |
|---|---|
| N no divisible por P (N=1003, P=4) | OK, partículas conservadas |
| Más procesos que partículas (N=3, P=4) ⇒ rangos vacíos | OK, sin crash |
| `--init twobody` distribuido (N se ajusta a 2) | OK, drift 3.2e-15 |
| Fugas de memoria (`leaks`, build secuencial) | 0 leaks for 0 total leaked bytes |
| Deadlocks | ninguno; solo se usan colectivos |

## Criterios de aceptación — Resumen

| Criterio | Umbral | Estado | Valor |
|---|---|---|---|
| `acc` con T hilos == con 1 hilo | bit a bit | **OK** | idéntico (test 9) |
| Σ n_local == N en todo paso | exacto | **OK** | test 10 |
| Checksum Σid, Σid² invariante | exacto | **OK** | test 11 |
| Splitters disjuntos y cubrientes | exacto | **OK** | test 12 |
| `acc` distribuido vs secuencial | < 1e-12 | **OK** | 0.0 exacto |
| Energía conservada (N=10K, 1000 pasos) | < 1% | **OK** | 6.51e-05 |
| Desbalance de partículas, P=8 | < 1.15 | **OK** | 1.0000 |
| Speedup OpenMP, 8 hilos | > 5× | **NO** | 3.94× (núcleos E) |
| Speedup MPI, 4 procesos | > 3× | **NO** | 2.80× (árbol duplicado) |
| N=500K con P=8 sin crash | — | **OK** | completa |
| Sin fugas ni deadlocks | limpio | **OK** | 0 leaks |
| Suite de validación | 13/13 | **OK** | 13/13 |

**10 de 12 criterios se cumplen.** Los dos que no se cumplen tienen causas identificadas y medidas, no son fallos de correctitud:

1. **OpenMP 3.94× < 5×** — causa de hardware. La fase paralela pura llega a 4.36×, lo que descarta Amdahl; el límite son los 4 núcleos *efficiency*, que rinden ~21% de un núcleo P. Con 4 hilos sobre núcleos P la eficiencia es 88%. Re-medir en pcunix.

2. **MPI 2.80× < 3×** — causa algorítmica y esperada. La fase de fuerzas sí cumple (3.48×); lo que arrastra el total es la construcción del árbol replicada en cada proceso. Es precisamente la limitación que el LET viene a eliminar en la semana 4.

## Argumentos CLI (nuevos)

```
--threads T            Hilos OpenMP (default: OMP_NUM_THREADS)
--rebalance-every K    Reparticion cada K pasos, solo MPI (default: 20)
--metrics              Desglose de tiempo por fase, solo MPI
```

Compilación y ejecución:

```bash
make                                  # produce nbody y nbody_mpi
./nbody --validate                    # 9 tests
mpirun -np 4 ./nbody_mpi --validate   # 13 tests
OMP_NUM_THREADS=2 mpirun -np 4 ./nbody_mpi -n 100000 -t 0.01 --metrics
```

En pcunix hay que cargar el módulo y forzar TCP antes (ver [mpi-setup.md](mpi-setup.md)):

```bash
module load mpi/mpich-x86_64
export FI_PROVIDER=tcp
```

## Pendiente de la semana

- **Corrida distribuida en varias máquinas** (`-hosts pcunix40,pcunix42`): no se pudo verificar en el portátil. Requiere las claves SSH aceptadas según mpi-setup.md.
- **Re-medir todos los benchmarks en pcunix**, con núcleos homogéneos y `OMP_PROC_BIND=close`, `OMP_PLACES=cores`.
- `--method direct` no está soportado en la versión distribuida (su rol es ser oráculo secuencial); el binario avisa y usa Barnes-Hut.

## Siguiente paso: Semana 4

Los datos de esta semana dejan el argumento para el LET cuantificado y, sobre todo, **corregido respecto de lo que anticipaba el plan**: el problema no es el volumen de comunicación (el `Allgatherv` es el 1.3% del tiempo) sino la **construcción del árbol replicada**, que crece con P y hace que agregar procesos empeore el rendimiento a partir de 4.

`replicate_particles` está aislada en una única función precisamente para que la semana 4 la sustituya por el intercambio de árboles localmente esenciales tocando un solo punto del código.

Segundo insumo para la semana 4: el desbalance de **trabajo** (1.11) es notoriamente mayor que el de **partículas** (1.00), así que el rebalanceo dinámico debería repartir por costo estimado del recorrido y no por conteo. La bisección por histograma de `domain_partition` ya sirve tal cual para eso: basta cambiar el peso de cada partícula de 1 a su costo medido.
