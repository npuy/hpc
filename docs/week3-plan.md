# Semana 3 — Plan Detallado: Versión Híbrida Básica (MPI + OpenMP)

## Objetivo

Pasar de la versión secuencial Barnes-Hut a una **versión híbrida MPI + OpenMP**:
las partículas se reparten entre procesos MPI por rangos de clave Morton, cada
proceso paraleliza su cálculo de fuerzas con OpenMP, y las partículas que salen
de su rango **migran** al proceso que corresponde.

Corresponde a las **etapas 5, 6 y 7** del plan general:

| Etapa | Objetivo       | Validación                    |
| ----- | -------------- | ----------------------------- |
| 5     | OpenMP         | Mismo resultado, speedup      |
| 6     | MPI estático   | Conservación de partículas    |
| 7     | Migración      | Ninguna partícula perdida     |

**Precondición (semanas 1-2, ya cumplidas):** Barnes-Hut secuencial validado,
8/8 tests PASS, speedup 14.8× vs directo en N=100K, sin fugas de memoria.

## Decisión clave de la semana: *replicación* en lugar de LET

El LET (Locally Essential Tree) es **etapa 8, semana 4**. Para tener una versión
híbrida funcionando esta semana sin depender de él, se usa un paso intermedio:

> Cada proceso **replica** el arreglo global de partículas vía `MPI_Allgatherv`,
> construye el **árbol global completo** en local, pero **solo calcula fuerzas
> para su propia porción** de partículas.

Ventajas:

- Las fuerzas son **numéricamente idénticas** a la versión secuencial ⇒ la
  validación es trivial (comparar contra el binario secuencial).
- El paralelismo de cómputo es real: el trabajo O(N log N) del solver se divide
  entre P procesos.
- Se ejercitan de verdad el particionamiento, la migración y el tipo MPI, que
  son exactamente las piezas que sobreviven a la semana 4.

Limitación conocida (y esperada): memoria y comunicación son O(N) por proceso, y
la construcción del árbol se repite en cada rango. Esto **acota la escalabilidad**
y es precisamente lo que el LET viene a resolver en la semana 4. Hay que
**instrumentarlo y reportarlo**, no esconderlo: el informe gana mucho mostrando
la curva de escalabilidad de la versión replicada y luego la mejora del LET.

## Estado de partida (qué existe)

```
src/
  vec3.h              -- macros vectoriales 3D              [reutilizar]
  particle.h / .c     -- struct Particle, init, I/O          [MODIFICAR: agregar id]
  force.h / .c        -- O(N^2) + energias                   [oraculo + OpenMP]
  leapfrog.h / .c     -- integrador KDK                      [OpenMP]
  morton.h / .c       -- morton_encode, morton_sort          [base del particionamiento]
  octree.h / .c       -- octree array-pool                   [reutilizar sin cambios]
  barnes_hut.h / .c   -- fuerza O(N log N)                   [OpenMP + variante por rango]
  validation.h / .c   -- 8 tests secuenciales                [extender: tests 9-13]
  main.c              -- loop secuencial                     [guardas #ifdef USE_MPI]
Makefile              -- targets nbody y nbody_mpi           [agregar -fopenmp y -DUSE_MPI]
```

## Estructura de archivos objetivo (nuevos)

```
src/
  mpi_types.h / .c    -- MPI_Datatype para Particle, init/free
  domain.h / .c       -- particionamiento por rangos Morton (biseccion por histograma)
  migration.h / .c    -- deteccion y reparto de particulas fuera de rango (Alltoallv)
  metrics.h / .c      -- timers por fase (compute / tree / mpi / migracion)
docs/specs/
  domain.md
  migration.md
  metrics.md
```

## Cambio de estructura: `Particle` necesita un `id`

Para poder verificar que **ninguna partícula se pierde ni se duplica** tras la
migración hace falta identidad estable. Se agrega un campo:

```c
typedef struct {
    double   pos[3];
    double   vel[3];
    double   acc[3];
    double   mass;
    uint64_t morton;
    int64_t  id;      /* NUEVO: identidad global estable, asignada en init */
} Particle;
```

- `init_plummer` / `init_uniform_cube` / `init_two_body` asignan `id = i`.
- El CSV suma una columna `id` (actualizar `write_particles_csv` y
  `read_particles_csv` de forma coherente).
- El checksum de validación es `Σ id` y `Σ id²` sobre todos los procesos: si
  ambos coinciden con los valores iniciales, no hubo pérdidas ni duplicados.

## Reparto de trabajo

- **Estudiante A → MPI y particiones:** `mpi_types.h/.c`, `domain.h/.c`,
  esqueleto MPI de `main.c`, `metrics.h/.c`.
- **Estudiante B → OpenMP y migración:** paralelización de `barnes_hut.c`,
  `force.c`, `leapfrog.c`, `morton.c`; `migration.h/.c`; tests 9-13.

Interfaz a acordar el **Día 1**: el struct `Domain` y las firmas de
`domain_partition`, `migrate_particles` y `compute_forces_bh_range`. Con eso
ambos pueden avanzar sin bloquearse.

---

## Día 1: OpenMP (rápido y de alto retorno)

### Estudiante B

El recorrido del árbol es **trivialmente paralelo**: cada partícula lee un árbol
inmutable y escribe solo su propio `acc[]`. No hay carreras ni acumulación
compartida.

```c
void compute_forces_bh(const Octree *t, Particle *p, int n,
                       double theta, double softening) {
    double eps2 = softening * softening, theta2 = theta * theta;
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; i++) {
        double acc[3] = {0, 0, 0};
        accumulate_force(t, t->root, i, p, theta2, eps2, acc);
        VEC3_COPY(p[i].acc, acc);
    }
}
```

**Por qué `schedule(dynamic, 64)`:** el costo por partícula no es uniforme — una
partícula en una región densa abre muchos más nodos que una en la periferia. Con
`static` los hilos que tocan zonas vacías terminan primero y quedan ociosos. El
chunk de 64 amortiza el costo de sincronización del scheduler.

Bucles a paralelizar:

| Función                       | Directiva                                        |
| ----------------------------- | ------------------------------------------------ |
| `compute_forces_bh`           | `parallel for schedule(dynamic, 64)`             |
| `potential_energy_bh`         | `parallel for schedule(dynamic, 64) reduction(+:U)` |
| `compute_forces_direct`       | `parallel for schedule(static)` (costo uniforme) |
| `kinetic_energy`              | `parallel for reduction(+:K)`                    |
| `leapfrog_kick` / `_drift`    | `parallel for schedule(static)`                  |
| cálculo de claves en `morton_sort` | `parallel for schedule(static)`             |

**No paralelizar:** `octree_build` (inserción con `realloc` sobre un pool
compartido — requiere árbol bottom-up por Morton, fuera de alcance) ni el `qsort`
de `morton_sort`. Ambos quedan como parte serial → **medirlos aparte**, porque son
el techo de Amdahl de esta etapa.

> ⚠ **`compute_forces_direct` NO es paralelizable tal cual — confirmado.**
> `src/force.c:28-30` aprovecha la 3ª ley de Newton acumulando en `p[j].acc`
> dentro del bucle interno. Con `parallel for` sobre `i`, dos hilos distintos
> escriben el mismo `p[j].acc` ⇒ **carrera de datos**, resultados no
> deterministas y silenciosamente incorrectos.
>
> Solución: agregar `compute_forces_direct_par` sin simetría (cada `i` recorre
> todos los `j != i` y escribe solo `p[i].acc`). Cuesta 2× flops pero es
> embarazosamente paralela y **bit a bit reproducible**. Conservar la versión
> simétrica intacta como oráculo secuencial: es la referencia de los tests 1-8 y
> no debe tocarse.

**Determinismo:** como cada iteración escribe solo `p[i].acc` y el orden de las
sumas *dentro* de una partícula no cambia, el resultado es **bit a bit idéntico**
al secuencial para cualquier número de hilos. Esto se testea (test 9). Las
reducciones de energía sí dependen del orden ⇒ tolerancia ~1e-12, no bit a bit.

### Nota de compilación local (macOS)

El `gcc` de Apple es clang y no trae OpenMP:

```bash
brew install libomp
# compilar con:
clang -Xpreprocessor -fopenmp -lomp ...
# o directamente con gcc real:
brew install gcc && make CC=gcc-14
```

En las máquinas de fing (`pcunix*`) `gcc -fopenmp` funciona sin ceremonia.

### Entregable Día 1 (B)

- `make` con `-fopenmp`; `OMP_NUM_THREADS=8 ./nbody -n 100000` corre.
- Speedup medido con 1/2/4/8 hilos y `acc` idéntico bit a bit al de 1 hilo.

---

## Día 2-3: MPI estático y particionamiento

### Estudiante A

#### 1. Tipo MPI para `Particle`

```c
void mpi_particle_type_init(MPI_Datatype *type);
void mpi_particle_type_free(MPI_Datatype *type);
```

Usar `MPI_Type_create_struct` con los bloques reales (`10×MPI_DOUBLE`,
`1×MPI_UINT64_T`, `1×MPI_INT64_T`) y `MPI_Type_create_resized` con
`sizeof(Particle)` como extent, para que el padding del compilador no rompa el
envío de arreglos. **No** usar `MPI_BYTE` sobre el struct crudo: funciona solo si
todos los nodos son homogéneos y oculta errores de alineación.

#### 2. Descomposición del dominio

```c
typedef struct {
    uint64_t min_key;    /* inicio del rango de claves (inclusive) */
    uint64_t max_key;    /* fin del rango (exclusive) */
    int      rank;
    int      nprocs;
    double   gmin[3], gmax[3];  /* bounding box global (para codificar Morton) */
} Domain;

/* Calcula el bounding box global con MPI_Allreduce(MIN/MAX). */
void domain_global_bounds(const Particle *p, int n_local, double gmin[3], double gmax[3]);

/* Reparte el espacio de claves entre nprocs de modo que cada rango
   contenga aproximadamente N/P particulas. Colectivo. */
void domain_partition(Domain *d, const Particle *p, int n_local, int64_t n_global);

/* Rango al que pertenece una clave (busqueda binaria sobre los splitters). */
int  domain_owner(const Domain *d, const uint64_t *splitters, uint64_t key);
```

**Algoritmo de particionamiento — bisección por histograma:**

Repartir el espacio de claves en partes iguales daría un desbalance enorme
(Plummer es fuertemente no uniforme: casi todas las partículas caen en unas pocas
celdas centrales). Hay que repartir por **conteo de partículas**, no por volumen:

```
para cada splitter objetivo k = 1..P-1:   # objetivo: k*N/P particulas a la izquierda
    lo = 0, hi = 2^63
    repetir ~63 veces (biseccion sobre el rango de claves):
        mid = (lo + hi) / 2
        c_local  = # de particulas locales con key < mid     # busqueda binaria, O(log n)
        MPI_Allreduce(&c_local, &c_global, 1, MPI_LONG_LONG, MPI_SUM, comm)
        si c_global < k*N/P: lo = mid
        si no:               hi = mid
    splitter[k] = lo
```

Optimización obligatoria: **bisecar los P-1 splitters simultáneamente** en un
único `MPI_Allreduce` de P-1 enteros por iteración, no P-1 Allreduce separados.
Así el costo total es ~63 colectivos, no 63·(P-1).

> **Por qué este algoritmo y no un sample sort:** es más simple de depurar, no
> necesita ordenar globalmente, y —lo importante— **es exactamente el mismo
> procedimiento que usará el rebalanceo dinámico de la semana 4**. Se escribe una
> vez y se reusa: en la semana 4 basta con volver a llamarlo cada K pasos.

#### 3. Esqueleto MPI en `main.c`

Todo el código MPI va bajo `#ifdef USE_MPI` para que `nbody` (secuencial) y
`nbody_mpi` compilen del mismo fuente, como ya prevé el Makefile.

```c
/* --- setup (una vez) --- */
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

if (rank == 0) generar N particulas e inicializar ids;
repartir N/P particulas a cada rango (Scatterv inicial);
domain_global_bounds(...);
morton_sort local;
domain_partition(...);          /* calcula splitters */
migrate_particles(...);         /* deja cada particula en su dueno */

/* --- bucle temporal --- */
for (cada paso) {
    leapfrog_kick(local, n_local, dt/2);
    leapfrog_drift(local, n_local, dt);

    recomputar claves Morton (las particulas se movieron);
    migrate_particles(...);                       /* etapa 7 */

    MPI_Allgatherv(local -> global, n_global);    /* replicacion (temporal) */
    Octree *tree = octree_build(global, n_global);
    octree_compute_mass(tree, global);
    compute_forces_bh_range(tree, global, n_global, off, off + n_local, theta, eps);
    copiar acc de global[off..off+n_local) a local;
    octree_free(tree);

    leapfrog_kick(local, n_local, dt/2);

    /* cada M pasos: energia global por Allreduce + chequeo de conservacion */
}
MPI_Finalize();
```

**Invariante que hace que esto funcione:** cada rango mantiene sus partículas
ordenadas por clave, y los rangos de claves son disjuntos y crecientes con el
rank. Entonces la concatenación que produce `Allgatherv` **ya está globalmente
ordenada por Morton** — sin ningún sort adicional. Vale la pena aseverarlo con un
`assert` en modo debug.

#### 4. Variante por rango en `barnes_hut.h`

```c
void compute_forces_bh_range(const Octree *t, Particle *p, int n,
                             int i0, int i1, double theta, double softening);
```

Idéntica a `compute_forces_bh` pero el bucle va de `i0` a `i1`. La versión
completa pasa a ser un wrapper con `(0, n)` — sin duplicar el código del
traversal. El pragma OpenMP vive en esta función: **ahí se compone el híbrido**
(MPI reparte `[i0,i1)` entre procesos, OpenMP reparte ese rango entre hilos).

### Entregable Día 3 (A)

- `mpirun -np 4 ./nbody_mpi -n 10000` corre y la suma de `n_local` es N.
- Los splitters producen un desbalance `max(n_local)/avg(n_local) < 1.1` en
  Plummer con P=4 y P=8.

---

## Día 3-4: Migración

### Estudiante B

```c
/* Reubica las particulas cuya clave salio del rango local.
   Devuelve el nuevo n_local; puede reasignar *p (realloc). */
int migrate_particles(Particle **p, int n_local, int *capacity,
                      const Domain *d, const uint64_t *splitters,
                      MPI_Datatype ptype, MPI_Comm comm);
```

Algoritmo (personalizado, colectivo, sin riesgo de deadlock):

```
1. recomputar clave Morton de cada particula local con el bbox GLOBAL
2. para cada particula: dest = domain_owner(splitters, key)
3. ordenar/particionar el arreglo local por dest  -> send_counts[P], send_displs[P]
4. MPI_Alltoall(send_counts -> recv_counts)
5. calcular recv_displs; asegurar capacidad del buffer de recepcion
6. MPI_Alltoallv(sendbuf -> recvbuf, ptype)
7. reemplazar el arreglo local por lo recibido; reordenar por clave Morton
8. n_local = suma(recv_counts)
```

Puntos finos:

- **El bounding box global cambia** cuando el sistema se expande o colapsa. Si se
  recodifica Morton con un bbox distinto al usado para calcular los splitters, la
  asignación es inconsistente. Decisión recomendada: **recalcular el bbox global
  y los splitters juntos**, cada K pasos (`--rebalance-every`, default 20), y
  entre medio usar el bbox congelado. Documentar la decisión.
- ⚠ **`morton_encode` no clampea — bug confirmado, hay que arreglarlo.** En
  `src/morton.c:24-26` la coordenada normalizada se castea a `uint64_t` sin
  acotarla. Con el bbox congelado, una partícula que se sale del box produce:
  - si `x > max[0]`: `ix > 2^21-1`, y el `v &= 0x1fffff` de `expand_bits`
    **envuelve** el valor ⇒ una partícula del extremo derecho recibe una clave
    del extremo izquierdo y migra al proceso equivocado;
  - si `x < min[0]`: el double es negativo y su conversión a `uint64_t` es
    **comportamiento indefinido** en C.

  En la semana 2 esto no se notaba porque el bbox se recalculaba cada paso y
  siempre contenía todo. Con bbox congelado entre reparticiones deja de ser
  cierto. Fix: acotar antes del cast (`if (fx < 0) fx = 0; if (fx > range) fx = range;`)
  sobre las tres coordenadas. Agregar un test que codifique un punto fuera del
  box y verifique que cae en la celda del borde.
- `MPI_Alltoallv` es colectivo y no puede deadlockear, a diferencia de un patrón
  `Send`/`Recv` armado a mano. Es la razón para preferirlo aquí.
- El caso `n_local == 0` en algún rango es legal y debe funcionar (buffers
  vacíos, `Allgatherv` con count 0). Testearlo explícitamente.

### Estudiante B — Extensión de la suite de validación

Los tests 9-13 se agregan a `validation.c`. Los que necesitan MPI van bajo
`#ifdef USE_MPI` y se corren con `mpirun -np 4 ./nbody_mpi --validate`.

**Test 9: Determinismo de OpenMP**
- N=5000, calcular `acc` con `OMP_NUM_THREADS=1` y luego con 2, 4, 8 hilos
  (usando `omp_set_num_threads`).
- Criterio: **igualdad bit a bit** (`memcmp` de los `acc[]`). Cualquier
  diferencia indica una carrera o una acumulación compartida.

**Test 10: Conservación de partículas**
- Tras 100 pasos con migración, `MPI_Allreduce(SUM, n_local) == N`.
- Criterio: igualdad exacta.

**Test 11: Identidad de partículas (sin pérdidas ni duplicados)**
- Checksum global `Σ id` y `Σ id²` antes y después de 100 pasos.
- Criterio: ambos idénticos. (`Σ id` sola no detecta un intercambio de dos ids;
  el segundo momento sí.)

**Test 12: Particiones disjuntas y cubrientes**
- Los splitters son estrictamente crecientes, `splitter[0] == 0`,
  `splitter[P] == UINT64_MAX`, y toda partícula local satisface
  `min_key <= key < max_key`.

**Test 13: Fuerzas MPI == fuerzas secuenciales**
- Misma semilla, mismo N, mismo θ; correr `./nbody` y `mpirun -np 4 ./nbody_mpi`.
- Comparar `acc` partícula por partícula (emparejando por `id`).
- Criterio: error relativo < 1e-12. Con el esquema de replicación el árbol es el
  mismo en todos los rangos, así que debería ser **bit a bit**; se deja 1e-12 de
  margen por si el orden de `Allgatherv` reordena algo.

---

## Día 5: Métricas, benchmarks y reporte

### Ambos

#### 1. Módulo de métricas

```c
typedef struct {
    double tree_time;       /* construccion del octree + centros de masa */
    double force_time;      /* traversal Barnes-Hut */
    double integrate_time;  /* leapfrog */
    double comm_time;       /* Allgatherv */
    double migrate_time;    /* Alltoallv + reordenamiento */
    double total_time;
    int64_t n_local;        /* para medir desbalance */
    int64_t steps;
} Metrics;

void metrics_reset(Metrics *m);
void metrics_report(const Metrics *m, MPI_Comm comm);  /* min/max/avg por rango */
```

Usar `MPI_Wtime()` (no `clock_gettime`) en el camino MPI. `metrics_report` hace
`Allreduce` con `MIN`, `MAX` y `SUM` de cada campo: el cociente **max/avg de
`force_time`** es la métrica de desbalance que justifica el rebalanceo de la
semana 4.

#### 2. Nuevos flags de CLI

| Flag                 | Descripción                                | Default |
| -------------------- | ------------------------------------------ | ------- |
| `--threads T`        | hilos OpenMP (si no, usa `OMP_NUM_THREADS`)| auto    |
| `--rebalance-every K`| repartición de splitters cada K pasos      | 20      |
| `--metrics`          | imprime el desglose por fase al final      | off     |

#### 3. Experimentos

Ejecutar en `pcunix*` (ver [docs/mpi-setup.md](mpi-setup.md): `module load
mpi/mpich-x86_64` y `export FI_PROVIDER=tcp`).

**a) Escalabilidad OpenMP (1 proceso, N=100K, 10 pasos):** 1, 2, 4, 8, 16 hilos.
Reportar speedup, eficiencia, y el **porcentaje de tiempo en la parte serial**
(construcción del árbol) — es el que predice dónde se aplana la curva.

**b) Escalabilidad MPI (1 hilo, N=100K, 10 pasos):** 1, 2, 4, 8 procesos.
Desglosar cómputo vs `Allgatherv`: se espera que la comunicación empiece a
dominar, y ese es justamente el argumento a favor del LET.

**c) Híbrido:** combinaciones a producto constante (ej. 16 unidades: 1×16, 2×8,
4×4, 8×2, 16×1) para encontrar el reparto óptimo procesos/hilos.

**d) Escalabilidad débil:** N/P = 25K constante (25K×1, 50K×2, 100K×4, 200K×8).
Con replicación se espera que **degrade** al crecer P — medirlo y explicarlo.

**e) Distribuido:** al menos una corrida en 2 máquinas con `-hosts`, para
verificar que funciona fuera de una sola máquina (requiere las claves SSH ya
aceptadas, ver mpi-setup.md).

---

## Criterios de aceptación

| Test / Criterio                                              | Umbral                    |
| ------------------------------------------------------------ | ------------------------- |
| `acc` con T hilos == `acc` con 1 hilo                        | bit a bit                 |
| Suma de `n_local` sobre los rangos == N                      | exacto, todo paso         |
| Checksum `Σ id` y `Σ id²` invariante tras 100 pasos          | exacto                    |
| Splitters disjuntos, crecientes y cubrientes                 | exacto                    |
| `acc` distribuido (P=4) vs secuencial                        | error rel. < 1e-12        |
| Energía conservada en la corrida híbrida (N=10K, 1000 pasos) | drift < 1%                |
| Desbalance `max(n_local)/avg(n_local)` con Plummer, P=8      | < 1.15                    |
| Speedup OpenMP con 8 hilos (N=100K)                          | > 5×                      |
| Speedup MPI con 4 procesos (N=100K, cómputo)                 | > 3×                      |
| `mpirun -np 8 ./nbody_mpi -n 500000` completa sin crash      | ---                       |
| Sin fugas de memoria ni deadlocks                            | `valgrind` / sin cuelgues |
| Suite completa                                               | 13/13 PASS                |

## Parámetros recomendados

| Parámetro          | Valor dev | Valor test |
| ------------------ | --------- | ---------- |
| N                  | 10K       | 100K–500K  |
| procesos MPI       | 2         | 4–8        |
| hilos OpenMP       | 2         | 4–8        |
| theta              | 0.5       | 0.5        |
| dt                 | 0.001     | 0.0005     |
| softening          | 0.01      | 0.01       |
| rebalance-every    | 20        | 20         |

## Notas técnicas

- **Un solo `main.c`:** guardas `#ifdef USE_MPI` en vez de un `main_mpi.c`
  paralelo. Evita que las dos versiones se desincronicen. OpenMP no necesita
  guardas explícitas (`#pragma omp` se ignora sin `-fopenmp`), pero si se llama a
  `omp_get_num_threads()` hay que protegerlo con `#ifdef _OPENMP`.
- **Makefile:** `nbody` con `-fopenmp`; `nbody_mpi` con `-fopenmp -DUSE_MPI` y
  `mpicc`. Agregar los cuatro `.c` nuevos a `SRC`.
- **`Allgatherv` es temporal por diseño.** Aislar la llamada detrás de una función
  (`exchange_particles(...)`) para que en la semana 4 se sustituya por el
  intercambio LET tocando un solo lugar.
- **Orden de las fases:** migrar **antes** del `Allgatherv`, no después. Si se
  migra al final del paso, el árbol de ese paso se construye con partículas que
  ya no pertenecen a nadie y el offset `[i0,i1)` deja de ser contiguo.
- **`MPI_Alltoallv` sobre `Send`/`Recv`:** colectivo, sin deadlocks, y el runtime
  elige el algoritmo (bruck / pairwise) según el tamaño. Un patrón manual es más
  código y más riesgo, sin beneficio a esta escala.
- **Afinidad de hilos:** en las máquinas de fing, `OMP_PROC_BIND=close` y
  `OMP_PLACES=cores` evitan que el SO migre hilos entre núcleos y ensucie las
  mediciones. Documentar la configuración usada en el reporte.
- **Sobresuscripción:** con P procesos × T hilos > núcleos físicos los tiempos se
  degradan de forma no lineal. Verificar `nproc` en la máquina antes de diseñar
  la grilla de experimentos híbridos.
- **Medir siempre con `MPI_Barrier` antes de arrancar el timer** de una fase
  colectiva; si no, el tiempo de comunicación absorbe el desbalance de cómputo y
  se malinterpreta el resultado.

## Riesgos de la semana

| Riesgo                                    | Mitigación                                                |
| ----------------------------------------- | --------------------------------------------------------- |
| Carrera en OpenMP (acc compartida)        | Test 9 bit a bit; revisar simetría en `compute_forces_direct` |
| Partículas perdidas en la migración       | Checksum `Σ id` / `Σ id²` cada paso en modo debug          |
| Splitters inconsistentes (bbox cambiante) | Congelar bbox entre reparticiones; clampear en `morton_encode` |
| Deadlock MPI                              | Solo colectivos (`Alltoallv`, `Allgatherv`, `Allreduce`)   |
| Desbalance por Plummer no uniforme        | Particionamiento por conteo (bisección), no por volumen    |
| Padding del struct en el tipo MPI         | `Type_create_struct` + `Type_create_resized`, nunca `MPI_BYTE` |
| `Allgatherv` domina y mata el speedup     | Instrumentar `comm_time`; es el argumento para el LET (sem. 4) |
| OpenMP no compila en macOS                | Desarrollar con `gcc-14` de brew o directo en `pcunix*`    |

## Entregables de la semana

1. `src/mpi_types.h/.c`, `src/domain.h/.c`, `src/migration.h/.c`,
   `src/metrics.h/.c` compilando sin warnings.
2. OpenMP en `barnes_hut.c`, `force.c`, `leapfrog.c`, `morton.c`.
3. `Particle` con `id`; CSV actualizado.
4. `main.c` híbrido con `#ifdef USE_MPI`; Makefile con `nbody` y `nbody_mpi`.
5. Tests 9-13 agregados (`--validate` pasa 13/13, secuencial y con `mpirun`).
6. `docs/specs/domain.md`, `migration.md`, `metrics.md`.
7. `docs/week3-report.md` con las tablas de escalabilidad (a-e) y el análisis del
   costo de `Allgatherv` que motiva el LET de la semana 4.
