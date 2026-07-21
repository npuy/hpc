# Semana 4 — Plan Detallado: LET, Balance Dinámico y Sistema Final

## Objetivo

Eliminar la replicación global de la semana 3 y sustituirla por el intercambio de
**árboles localmente esenciales (LET)**, agregar **rebalanceo dinámico por costo**
y producir los **experimentos y el informe final**.

Corresponde a las **etapas 8 y 9** del plan general:

| Etapa | Objetivo   | Validación                      |
| ----- | ---------- | ------------------------------- |
| 8     | LET        | Comparación contra árbol global |
| 9     | Rebalanceo | Reducción del desbalance        |

## Punto de partida: lo que midió la semana 3

Dos números del [reporte de la semana 3](week3-report.md) dirigen todo el diseño
de esta semana, y **el primero corrige lo que anticipaba el plan**:

1. **El cuello de botella no es la comunicación.** El `Allgatherv` es el **1.3%**
   del tiempo total con P=8. Lo caro es la **construcción del árbol replicada**:
   `tree_time` crece 4.8× de P=1 a P=8 (0.235 → 1.133 s) porque cada proceso
   construye el árbol sobre las N partículas. Por eso el speedup MPI se da vuelta
   después de P=4 (2.80× con 4, 2.61× con 8) y la escalabilidad débil cae al 24%.

   ⇒ **El LET no se justifica por ahorro de red sino por ahorro de cómputo
   redundante.** El objetivo medible de la semana es que `tree_time` deje de
   crecer con P.

2. **El desbalance de trabajo (1.11) es mayor que el de partículas (1.00).**
   `domain_partition` reparte los conteos casi perfectamente, pero una partícula
   en zona densa abre más nodos que una en la periferia.

   ⇒ **El rebalanceo debe repartir por costo estimado, no por conteo.** La
   bisección por histograma sirve tal cual: solo cambia el peso de cada partícula
   de 1 a su costo medido.

## Estado de partida (qué existe y qué hay que tocar)

```
src/
  octree.h / .c       -- pool de nodos, box calculado de las particulas  [MODIFICAR: box explicito]
  barnes_hut.h / .c   -- traversal + OpenMP                              [MODIFICAR: contador de costo]
  domain.h / .c       -- biseccion por histograma de CONTEOS             [MODIFICAR: pesos]
  migration.h / .c    -- Alltoallv + replicate_particles                 [replicate -> queda como opcion]
  metrics.h / .c      -- timers por fase                                 [MODIFICAR: fases LET]
  mpi_types.h / .c    -- tipo derivado de Particle                       [MODIFICAR si se agrega campo]
  particle.h / .c     -- struct con id                                   [MODIFICAR: campo work]
  main.c              -- run_mpi con replicacion                         [MODIFICAR: --exchange]
  validation.h / .c   -- 13 tests                                        [extender: tests 14-18]
```

Nuevos:

```
src/
  let.h / .c          -- seleccion, empaquetado e importacion del arbol localmente esencial
docs/specs/
  let.md
  balance.md
```

## Reparto de trabajo

- **Estudiante A → LET:** `let.h/.c`, `octree_build_box`, integración en `main.c`,
  tests 14-16.
- **Estudiante B → Balance:** contador de costo en `barnes_hut.c`, campo `work`,
  `domain_partition_weighted`, métricas nuevas, tests 17-18.
- **Ambos → Día 5:** experimentos e informe final.

Interfaz a acordar el **Día 1**: firma de `let_exchange` y del contador de costo.
Con eso ambos avanzan sin bloquearse.

---

## Día 1: Preparación — el árbol necesita una caja compartida

### Estudiante A

**El problema:** `octree_build` llama internamente a `octree_bounds`, que calcula
el cubo a partir de las partículas recibidas. Con replicación todos los procesos
pasaban el mismo arreglo y obtenían el mismo cubo. Con LET cada proceso construye
sobre sus N/P locales ⇒ **cada proceso obtendría un cubo distinto**, las celdas no
alinearían entre procesos y las aproximaciones importadas no serían comparables.

Corrección mínima:

```c
/* Construye el arbol usando un cubo raiz EXPLICITO en vez de derivarlo de las
   particulas. Todos los procesos deben pasar el mismo cubo (d->gmin/gmax) para
   que sus celdas coincidan: es la precondicion del LET. */
Octree *octree_build_box(const Particle *p, int n,
                         const double min[3], const double max[3]);
```

`octree_build` pasa a ser el wrapper que calcula `octree_bounds` y delega. Cero
cambios en el resto del código secuencial.

> **Precondición nueva:** con caja fija, una partícula fuera del cubo global no
> tiene octante válido. Es el mismo escenario que ya obligó a acotar en
> `morton_encode` (bug 1 de la semana 3). `octree_build_box` debe **acotar la
> posición al insertar** (no modificar la partícula, solo el octante elegido) o
> bien afirmar la precondición con un `assert`. Decidir y documentar.

### Estudiante B

Instrumentar el costo del recorrido. `accumulate_force` pasa a llevar un contador
de interacciones (hojas visitadas + nodos aceptados):

```c
static void accumulate_force(..., double acc[3], int64_t *work);
```

y `compute_forces_bh_range` lo guarda en `p[i].work`. Es un incremento por
interacción: costo despreciable frente al `sqrt` que ya hay en cada una.

Campo nuevo en `Particle`:

```c
typedef struct {
    double   pos[3], vel[3], acc[3];
    double   mass;
    uint64_t morton;
    int64_t  id;
    int64_t  work;   /* NUEVO: interacciones del ultimo paso (peso del balance) */
} Particle;
```

⚠ Agregar el campo **obliga a actualizar `mpi_particle_type_init`** (un bloque más
de `MPI_INT64_T`). Si no se actualiza, la migración copia basura silenciosamente.
Es el riesgo más fácil de olvidar de la semana.

> **Alternativa considerada y descartada:** mantener los pesos en un arreglo
> paralelo. Obliga a migrarlo por separado y a mantener su permutación sincronizada
> con la de las partículas. Meter el peso adentro del struct lo hace viajar solo,
> a costa de 8 bytes por partícula.

### Entregable Día 1

- `octree_build_box` con `octree_build` como wrapper; suite 13/13 sin regresión.
- `p[i].work` poblado; tipo MPI actualizado; test 11 (checksum) sigue PASS.

---

## Día 2-3: LET (Locally Essential Tree)

### Estudiante A

#### Idea

El árbol localmente esencial de un proceso `r` es el subconjunto del árbol global
que `r` necesita para calcular las fuerzas de **sus** partículas con la misma
precisión que tendría con el árbol completo. Todo lo demás está lo bastante lejos
como para resumirse en unos pocos nodos agregados.

Cada proceso recorre **su** árbol local una vez por destino y decide, nodo por
nodo, si el destino puede quedarse con el resumen (masa + centro de masa) o
necesita el detalle.

#### 1. Cajas de dominio

El criterio de exportación necesita saber *dónde está* el destino. Un tramo Morton
no es una caja, así que cada proceso publica el AABB de sus partículas locales:

```c
/* Allgather de 6 doubles por proceso: P cajas. Despreciable frente a Alltoallv. */
void let_gather_domain_boxes(const Particle *local, int n_local,
                             double *boxes /* 6*P */, MPI_Comm comm);
```

Se recalcula **cada paso**, después de migrar y antes de exportar: las partículas
se movieron y una caja obsoleta puede dejar afuera interacciones necesarias.

#### 2. Criterio de exportación (MAC conservador)

Para el nodo `c` (lado `s`, centro de masa `cm`) y la caja destino `B_r`:

```
d_min = distancia minima de cm a la caja B_r   (0 si cm cae dentro de B_r)
aceptar si  s² < theta² · (d_min² + eps²)
```

**Por qué esto es correcto.** El receptor, al calcular la fuerza sobre su partícula
`i ∈ B_r`, aplicaría `s² < theta²·(d_i² + eps²)` con `d_i = |cm - pos_i| ≥ d_min`.
Como usamos la distancia *mínima* sobre toda la caja, si nosotros aceptamos el
resumen, el receptor también lo habría aceptado. **Nunca exportamos un resumen
donde el receptor quería detalle.** El caso inverso (abrimos un nodo que el
receptor habría aceptado) manda datos de más: cuesta ancho de banda, no precisión.

#### 3. Selección y empaquetado

```c
para cada destino r != rank:
    recorrer el arbol local desde la raiz:
        si nodo vacio (mass == 0):        podar
        si nodo interno y MAC(c, B_r):    exportar PSEUDO-PARTICULA (mass, cm)
        si nodo interno y no MAC:         descender a los 8 hijos
        si hoja:                          exportar las particulas REALES de la hoja
```

Las entradas exportadas se empaquetan como `Particle` con `pos = cm`,
`mass = node->mass`, `vel = acc = 0`, `id = -1` (marca de fantasma). Así se
reutiliza `ptype` y `MPI_Alltoallv` tal cual, sin un tipo MPI nuevo.

```c
/* Selecciona, intercambia e importa el LET. Los fantasmas quedan APPENDEADOS
   despues de las n_local locales en *p (que puede reallocarse). Retorna n_ghost. */
int let_exchange(Particle **p, int n_local, int *capacity,
                 const Octree *local_tree, const Domain *d,
                 double theta, double softening,
                 MPI_Datatype ptype, MPI_Comm comm);
```

#### 4. Cálculo de fuerzas con fantasmas

Los fantasmas viven en `[n_local, n_local + n_ghost)`. Las fuerzas se calculan en
`[0, n_local)` sobre un árbol que incluye a ambos:

```c
Octree *t = octree_build_box(p, n_local + n_ghost, d->gmin, d->gmax);
octree_compute_mass(t, p);
compute_forces_bh_range(t, p, n_local + n_ghost, 0, n_local, theta, softening);
```

Es exactamente la firma que ya existe: `compute_forces_bh_range` fue escrita en la
semana 3 para esto.

> **Optimización para el Día 4, no para el Día 2:** se construyen dos árboles por
> paso (el local para seleccionar, el fusionado para calcular). Como la caja raíz
> es la misma, se pueden **insertar los fantasmas en el árbol local ya construido**
> (`octree_insert_particles(t, p, n_local, n_local+n_ghost)`, ~15 líneas, hay que
> reallocar `t->next`). Ahorra una construcción completa. Hacerlo **después** de
> que la versión con dos árboles valide, no antes.

#### 5. Cuidados

- **`particle_checksum` y las energías deben ignorar los fantasmas.** Contar el
  rango `[0, n_local)` y nada más. Un fantasma en el checksum lo rompe con
  `id = -1`; un fantasma en la energía cinética la infla.
- **Ningún proceso se importa a sí mismo:** solo se exporta a `r != rank`, así que
  no hay riesgo de contar dos veces la propia masa.
- **`n_ghost` es dinámico y crece cuando el sistema se concentra.** Manejar la
  capacidad del buffer con la misma lógica de realloc que `migrate_particles`.
- **La suma de masas es un invariante barato:** `Σ mass` sobre locales + fantasmas
  de un proceso debe dar ≈ la masa total del sistema (el LET resume, no descarta).
  Es la primera verificación a escribir, y detecta casi cualquier bug de selección.

### ⚠ La validación cambia de naturaleza: se pierde la exactitud bit a bit

El test 13 de la semana 3 daba **0.0 exacto** porque todos los procesos construían
el mismo árbol global. **Con LET eso deja de ser cierto y es correcto que así sea**,
por una razón concreta:

> Un fantasma llega como punto de masa en `cm`. Al reinsertarlo en el árbol
> fusionado, puede quedar agrupado con otros fantasmas bajo un ancestro común que
> a su vez satisface el MAC ⇒ se aproxima **dos veces**. El error sigue acotado por
> θ pero ya no es cero.

Hay que reemplazar el criterio, no relajarlo a ojo. La escalera de validación:

| Nivel | Qué compara | Criterio |
|---|---|---|
| **θ = 0** | LET vs directo O(N²) | **bit a bit** — con θ=0 nada se aproxima, el LET degenera en replicación completa y el árbol fusionado tiene la misma estructura y el mismo orden de recorrido |
| **Referencia sin reagrupamiento** | árbol fusionado vs (traversal local + suma directa sobre la lista de fantasmas) | aísla y cuantifica *solo* el error de reagrupamiento |
| **Producción** | LET (θ=0.5) vs replicado (θ=0.5) | error relativo del orden del error propio de Barnes-Hut, no mayor |

El nivel θ=0 es el test fuerte: es exacto, es barato a N=2000, y detecta cualquier
interacción faltante. Escribirlo **primero**.

### Entregable Día 3 (A)

- `mpirun -np 4 ./nbody_mpi --exchange let -n 10000` corre; masa conservada.
- Test θ=0 bit a bit contra el directo, PASS.
- `n_ghost` reportado y **notoriamente menor que N**.

---

## Día 3-4: Balance dinámico por costo

### Estudiante B

La bisección por histograma de `domain_partition` ya hace lo correcto; solo hay
que cambiar qué cuenta. En vez de "cuántas partículas hay debajo de esta clave",
"cuánto **trabajo** hay debajo de esta clave":

```c
/* Igual que domain_partition pero reparte el trabajo (Σ p[i].work) en vez del
   conteo. Con weights=NULL o en el primer paso (work sin medir) cae al conteo. */
void domain_partition_weighted(Domain *d, const Particle *p, int n_local,
                               int64_t w_global, MPI_Comm comm);
```

Implementación: el arreglo local ya está ordenado por clave, así que se precomputa
una vez por llamada el **prefijo de `work`**; el conteo por candidato pasa de
"búsqueda binaria + índice" a "búsqueda binaria + lectura del prefijo". Mismo
O(log n) por candidato, misma cantidad de colectivos (~63), mismo código de
bisección simultánea.

Detalles:

- **El peso está retrasado un paso** (es el costo medido en el paso anterior). Es
  la práctica estándar y es válida porque la distribución cambia poco entre pasos
  consecutivos; conviene decirlo explícitamente en el informe.
- **Primer paso:** `work` vale 0 para todas. Detectar `w_global == 0` y caer al
  particionamiento por conteo.
- **Flag `--balance count|work`** (default `work`). Sin el flag no hay forma de
  medir la mejora, y la comparación de las dos curvas es uno de los resultados del
  informe.
- **Métrica de éxito:** el desbalance de `force_time` (max/avg), no el de conteo.
  Con pesos se espera que **el desbalance de partículas empeore** (los procesos con
  zonas densas reciben menos partículas) mientras el de trabajo mejora. Es la señal
  de que funciona, no un problema.

### Estudiante B — Tests 14-18

**Test 14: LET con θ=0 == directo (bit a bit)**
- N=2000, P=4, θ=0. Comparar `acc` partícula por partícula emparejando por `id`.
- Criterio: **igualdad bit a bit**. Es el test que garantiza que no falta ninguna
  interacción.

**Test 15: LET vs replicado a θ=0.5**
- Mismo N, semilla y θ; una corrida con `--exchange replicate` y otra con `let`.
- Criterio: error relativo máximo de `acc` **< 1e-3**, y el error del LET contra el
  directo O(N²) no mayor a **1.5×** el error del Barnes-Hut secuencial contra el
  directo. Reportar los tres números, no solo el veredicto.

**Test 16: Volumen del LET**
- Con N=100K y P=4/8, verificar `n_ghost < N/4` y registrar la relación
  `n_ghost / n_local` en función de P.
- Criterio: `n_ghost` **no crece linealmente con N** al mantener N/P fijo (es la
  propiedad que rescata la escalabilidad débil).

**Test 17: El rebalanceo por costo reduce el desbalance de trabajo**
- N=100K, P=8, 100 pasos; medir max/avg de `force_time` con `--balance count` y con
  `--balance work`.
- Criterio: desbalance de trabajo **< 1.05** con pesos (línea base medida en la
  semana 3: **1.1149** con conteo).

**Test 18: Conservación de energía con LET**
- N=10K, P=4×2 hilos, 1000 pasos, `--exchange let`.
- Criterio: drift **< 1%** (línea base con replicación: 6.51e-05). Un LET con
  interacciones faltantes se delata acá aunque el test 15 pase por poco.

---

## Día 5: Experimentos e informe final

### Ambos

Todo en `pcunix*` (ver [mpi-setup.md](mpi-setup.md): `module load
mpi/mpich-x86_64`, `export FI_PROVIDER=tcp`, `OMP_PROC_BIND=close`,
`OMP_PLACES=cores`).

> **Deuda de la semana 3:** los benchmarks se corrieron en un portátil Apple Silicon
> con 4 núcleos P + 4 E, lo que distorsiona todo lo que pase de 4 hilos. **Hay que
> re-medir la semana 3 completa en pcunix**, no solo la semana 4: sin línea base
> homogénea las tablas comparativas del informe no se sostienen. Presupuestar tiempo
> para esto el Día 5 temprano, no al final.

| # | Experimento | Qué debe mostrar |
|---|---|---|
| a | **Replicación vs LET**, N=100K, P=1..16, desglose por fase | `tree_time` deja de crecer con P; es **el** gráfico del informe |
| b | **Escalabilidad fuerte**, N=100K y 500K, P=1,2,4,8,16 | speedup > 3× con P=4 (la semana 3 se quedó en 2.80×) |
| c | **Escalabilidad débil**, N/P = 25K | eficiencia > 60% con P=8 (la semana 3 cayó a 24%) |
| d | **Híbrido**, P×T constante | el óptimo debería **dejar de ser 1×8**: con LET los procesos ya no duplican el árbol |
| e | **Balance**, `count` vs `work`, P=8 y 16 | desbalance de trabajo y tiempo total |
| f | **Volumen del LET**, n_ghost vs P y vs N | validación empírica de la escalabilidad |
| g | **Precisión**, θ ∈ {0.2, 0.5, 0.8, 1.0} | error vs tiempo, LET y secuencial superpuestos |
| h | **Distribuido**, `-hosts pcunix40,pcunix42` | deuda pendiente de la semana 3; al menos una corrida |

El experimento **(d)** es el que cierra el argumento: la conclusión de la semana 3
fue "conviene minimizar procesos MPI y maximizar hilos, porque los hilos comparten
el árbol y los procesos lo duplican". **Si el LET funciona, esa conclusión debe
invertirse.** Es la verificación más honesta de que la semana sirvió para algo.

### Informe final (`docs/final-report.md`)

Estructura sugerida: problema y algoritmo → arquitectura híbrida → descomposición
Morton → LET → balance dinámico → validación (18 tests) → resultados
experimentales → análisis crítico → conclusiones y trabajo futuro.

Dos cosas que el informe debe decir explícitamente, porque son el contenido real
del trabajo:

1. **El diagnóstico de la semana 3 refutó la hipótesis del plan** (se esperaba que
   dominara la comunicación; dominaba el cómputo redundante). Mostrar la medición
   que lo refuta.
2. **Los dos bugs de la semana 3** (`morton_encode` sin acotar, carrera en
   `compute_forces_direct`) y cómo los detectó la validación, no la inspección.

---

## Criterios de aceptación

| Criterio | Umbral |
|---|---|
| LET con θ=0 == directo O(N²) | bit a bit |
| LET (θ=0.5) vs replicado | error rel. < 1e-3 |
| Error LET vs directo / error BH secuencial vs directo | < 1.5× |
| Σ masa (locales + fantasmas) == masa total | < 1e-12 rel. |
| Checksum `Σ id`, `Σ id²` invariante (fantasmas excluidos) | exacto |
| Energía conservada, N=10K, 1000 pasos, LET | drift < 1% |
| `tree_time` con P=8 vs P=1 | **no crece** (semana 3: 4.8×) |
| Speedup MPI, P=4, N=100K | > 3× |
| Eficiencia débil, N/P=25K, P=8 | > 60% (semana 3: 24%) |
| Desbalance de trabajo con `--balance work`, P=8 | < 1.05 (semana 3: 1.115) |
| `n_ghost` con N=100K, P=8 | < N/4 |
| N=1M con P=16 sin crash | — |
| Corrida en ≥2 máquinas con `-hosts` | funciona |
| Sin fugas ni deadlocks | limpio |
| Suite completa | 18/18 PASS |

## Riesgos de la semana

| Riesgo | Impacto | Mitigación |
|---|---|---|
| **Fantasmas faltantes** (fuerzas mal por defecto, silenciosamente) | Alto | Test θ=0 bit a bit **primero**; invariante Σ masa cada paso en debug |
| Cajas raíz distintas entre procesos | Alto | `octree_build_box` con `d->gmin/gmax`; afirmarlo con `assert` |
| Tipo MPI desactualizado tras agregar `work` | Alto | Test 11 lo detecta; correrlo apenas se cambia el struct |
| Fantasmas contaminan checksum / energía / migración | Medio | Rangos explícitos `[0,n_local)`; nunca migrar `id == -1` |
| `n_ghost` explota al concentrarse el sistema | Medio | Reportarlo por paso; si crece, subir la frecuencia de repartición |
| Error de reagrupamiento mayor al esperado | Medio | Referencia sin reagrupamiento (dos pasadas) para separarlo del resto |
| Rebalanceo oscila (pesos retrasados) | Bajo | `--rebalance-every` ≥ 10; medir el costo de la repartición aparte |
| No queda tiempo para el informe | **Alto** | Congelar código el Día 4; el Día 5 es solo medir y escribir |

## Prioridades si falta tiempo

Del plan general: LET antes que balance. Concretamente:

1. LET correcto (tests 14 y 18) — sin esto la semana no tiene resultado.
2. Experimento (a): replicación vs LET, desglose por fase.
3. Escalabilidad fuerte y débil (b, c).
4. Balance por costo (test 17, experimento e).
5. Optimización de un solo árbol por paso.
6. Experimentos g, h.

**El balance dinámico es sacrificable; el LET no.** Si el Día 4 el LET todavía no
valida, abandonar el balance y dedicar los dos estudiantes al LET: un informe con
LET y sin balance es un trabajo completo, uno con balance y sin LET no.

## Parámetros recomendados

| Parámetro | Valor dev | Valor test |
|---|---|---|
| N | 10K | 100K–1M |
| procesos MPI | 4 | 8–16 |
| hilos OpenMP | 2 | 4–8 |
| theta | 0.5 | 0.5 |
| dt | 0.001 | 0.0005 |
| softening | 0.01 | 0.01 |
| rebalance-every | 20 | 10–20 |

## Argumentos CLI nuevos

```
--exchange replicate|let   Esquema de intercambio (default: let)
--balance count|work       Criterio de reparticion (default: work)
--theta-let X              Theta del criterio de exportacion (default: = --theta)
```

`--exchange replicate` **debe seguir funcionando**: es el oráculo de los tests 15-16
y la línea base de todos los gráficos comparativos. No borrar `replicate_particles`.

## Entregables de la semana

1. `src/let.h/.c`, `octree_build_box`, `domain_partition_weighted` sin warnings.
2. Campo `work` en `Particle` y tipo MPI actualizado.
3. `main.c` con `--exchange` y `--balance`; ambos caminos funcionando.
4. Tests 14-18 (`mpirun -np 4 ./nbody_mpi --validate` → 18/18).
5. `docs/specs/let.md`, `docs/specs/balance.md`.
6. `docs/week4-report.md` con los experimentos a-h.
7. `docs/final-report.md`: informe final del proyecto.
