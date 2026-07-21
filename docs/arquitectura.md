# Arquitectura: cómo funciona todo

Vista consolidada del simulador después de 5 semanas. Para el detalle de cada
pieza están las specs en `docs/specs/`; acá está cómo encajan.

## Qué hace el programa

Simula la evolución gravitatoria de N cuerpos en 3D. El costo está en calcular
las fuerzas: la suma directa es O(N²), y Barnes-Hut la baja a O(N log N)
aproximando grupos lejanos por su centro de masa. El paralelismo es híbrido:
**MPI** reparte partículas entre procesos, **OpenMP** reparte el trabajo de cada
proceso entre hilos.

## Los dos binarios

| Binario | Compilación | Rol |
|---|---|---|
| `nbody` | `gcc -fopenmp` | Secuencial + OpenMP. Es el **oráculo**: los tests comparan contra él. |
| `nbody_mpi` | `mpicc -fopenmp -DUSE_MPI` | Híbrido MPI + OpenMP. |

Un solo `main.c` con guardas `#ifdef USE_MPI`, para que las dos versiones no se
desincronicen.

## El paso de simulación

Integrador leapfrog en forma KDK (kick-drift-kick), que es simpléctico: conserva
la energía a largo plazo en vez de derivar.

```
leapfrog_kick(dt/2)            v += a·dt/2
leapfrog_drift(dt)             x += v·dt

[cada K pasos] REPARTICIÓN
    domain_global_bounds       bbox global por Allreduce MIN/MAX, congelado
    morton_sort                orden local (solo localidad de caché)
    domain_partition_orb       bisección recursiva del trabajo -> caja por proceso

MIGRACIÓN
    migrate_particles          dueño = domain_owner_pos(x), Alltoallv

ÁRBOL + LET
    octree_build_box           árbol local, caja raíz COMPARTIDA (gmin/gmax)
    octree_compute_mass        masas y centros de masa
    let_exchange               selección con MAC conservador + Alltoallv
                               -> fantasmas appendeados en [n_local, n_local+n_ghost)
    octree_insert_particles    fantasmas al árbol ya construido (no se reconstruye)
    octree_compute_mass        las masas cambiaron

FUERZAS
    compute_forces_bh_range(..., 0, n_local)     recorrido MAC, OpenMP

leapfrog_kick(dt/2)            v += a·dt/2
```

**El orden importa.** Migrar antes de construir el árbol, no después: si se
migrara al final, el árbol del paso se construiría con partículas que ya no
pertenecen a nadie.

## Las tres ideas que sostienen el diseño

### 1. Barnes-Hut: el criterio de apertura

Para cada partícula se recorre el árbol desde la raíz. Si una celda de lado `s`
está a distancia `d` y cumple `s/d < θ`, todo su contenido se aproxima por un
único cuerpo en su centro de masa; si no, se desciende a los 8 hijos.

Con `θ = 0` no se aproxima nada y el resultado coincide con el método directo.
Esa propiedad es la que hace posible el test más fuerte de la suite.

### 2. Descomposición del dominio: ORB

Cada proceso es dueño de una **caja** del espacio. Se obtiene bisecando
recursivamente el grupo de procesos: en cada nivel se corta por la dimensión más
larga, buscando el plano que reparte el **trabajo** por la mitad.

El árbol ORB tiene `2P-1` nodos y está replicado en todos los procesos — son
cientos de bytes, así que cada uno conoce la caja de todos los demás sin
comunicación.

> **Por qué cajas y no rangos de curva Morton.** Las semanas 3 y 4 usaron tramos
> contiguos de clave Morton. Funciona para migrar, pero **un tramo de curva no es
> una región compacta**: se midió que el AABB del dominio de un proceso cubría el
> 44,6% del espacio. Con dominios así el LET no puede resumir nada. Ver
> [specs/orb.md](specs/orb.md).

### 3. LET: cada proceso recibe solo lo que necesita

En vez de replicar las N partículas, cada proceso importa el subconjunto del árbol
global que necesita para calcular **sus** fuerzas con la misma precisión.

El criterio de exportación es un MAC **conservador**: para el nodo `c` y la caja
destino `B_r`, se usa la distancia **mínima** de `cm` a `B_r`. Como toda partícula
del destino está al menos a esa distancia, si nosotros aceptamos el resumen el
receptor también lo habría aceptado. Nunca se exporta un resumen donde el receptor
quería detalle.

Lo exportado viaja como `Particle`: los nodos aceptados como pseudo-partículas
(`pos = cm`, `id = -1`), las hojas como partículas reales. Así se reutiliza el
tipo MPI derivado sin empaquetado a mano.

**Los fantasmas son transitorios**: el `migrate_particles` del paso siguiente
reemplaza el arreglo local y los descarta solo. Todo lo que sea checksum, energía
o migración usa el rango `[0, n_local)`.

## Balance de carga

El desbalance de **partículas** y el de **trabajo** no son lo mismo: una partícula
en zona densa abre muchos más nodos que una en la periferia. La semana 3 midió
1.00 de desbalance de partículas contra 1.11 de trabajo.

`accumulate_force` cuenta las interacciones que evalúa y las guarda en
`Particle::work`. La bisección de ORB reparte esa magnitud en vez del conteo. El
peso vive dentro del struct (no en un arreglo paralelo) para que viaje solo en la
migración.

**Señal esperada y contraintuitiva:** al activar el reparto por costo, el
desbalance de *partículas* empeora mientras el de *trabajo* mejora. Es el punto,
no un problema.

## Mapa de módulos

```
src/
  vec3.h            macros vectoriales 3D
  particle.h/.c     struct Particle, generadores (plummer/uniform/twobody), CSV
  leapfrog.h/.c     integrador KDK
  force.h/.c        O(N²) directo (oráculo) + variante sin simetría paralelizable
  morton.h/.c       claves Z-order 21 bits/eje; hoy solo localidad de caché
  octree.h/.c       pool de nodos; build, build_box, insert_particles, compute_mass
  barnes_hut.h/.c   recorrido MAC, variante por rango, contador de costo, OpenMP
  domain.h/.c       ORB (árbol de cajas) + camino Morton legado
  migration.h/.c    Alltoallv por dueño; replicate (legado/oráculo); checksums
  let.h/.c          cajas de dominio, selección, intercambio, verificación de masa
  mpi_types.h/.c    tipo MPI derivado de Particle (4 bloques)
  metrics.h/.c      timers por fase, desbalances, volumen del LET
  validation.h/.c   19 tests
  main.c            CLI, run_sequential, run_mpi, build_step_tree
```

### `Particle`

```c
typedef struct {
    double   pos[3], vel[3], acc[3];
    double   mass;
    uint64_t morton;   /* clave Z-order */
    int64_t  id;       /* identidad estable; -1 marca un FANTASMA del LET */
    int64_t  work;     /* interacciones del último paso: peso del rebalanceo */
} Particle;
```

`id` sobrevive al reordenamiento y a la migración: es lo que permite verificar que
no se pierde ni se duplica ninguna partícula, y emparejar contra el oráculo
secuencial después de que la migración permute el arreglo.

> ⚠ Agregar un campo obliga a actualizar `mpi_particle_type_init`. Si no, la
> migración copia basura en silencio. El test 11 lo detecta.

## Caminos alternativos conservados

Cada esquema viejo sigue vivo detrás de un flag. **No es deuda técnica: es el
instrumento de medición.** Sin poder correr los dos caminos no hay forma de medir
la mejora, y todas las tablas comparativas de los reportes salen de ahí.

| Flag | Default | Alternativa |
|---|---|---|
| `--decomp` | `orb` | `morton` — rangos de clave (semanas 3-4) |
| `--exchange` | `let` | `replicate` — réplica global vía Allgatherv (semana 3) |
| `--balance` | `work` | `count` — reparto por conteo de partículas |

## Validación: 19 tests

La suite está construida en capas, de propiedades físicas a comparaciones contra
oráculos.

| # | Qué verifica | Criterio |
|---|---|---|
| 1-3 | Órbita de dos cuerpos, energía, momento | analítico / conservación |
| 4 | Codificación y orden Morton | 8 esquinas del cubo |
| 5-6 | Masa y centro de masa del árbol | vs suma directa |
| 7-8 | Fuerza Barnes-Hut y convergencia en θ | vs O(N²) |
| 9 | Determinismo de OpenMP | **bit a bit** con 1/2/4/8 hilos |
| 10-11 | Conservación de partículas e identidad | Σid y Σid² invariantes |
| 12 | Descomposición válida (ORB **y** Morton) | disjunta, cubriente, coherente |
| 13 | Fuerzas replicadas vs secuenciales | error < 1e-12 |
| 14 | **LET con θ=0 vs BH secuencial** | **bit a bit** |
| 15 | Precisión del LET a θ=0.5 | < 1e-3, y < 1.5× el error del BH |
| 16 | Crecimiento de `n_ghost` con N | exponente k < 0.9 |
| 17 | Rebalanceo por costo | desbalance de trabajo < 1.05 |
| 18 | Energía con LET, 200 pasos | drift < 1% |
| 19 | Volumen del LET, ORB vs Morton | ORB < 0.6× Morton |

**El test 14 es el guardián.** Con θ=0 el LET degenera en réplica completa, el
árbol fusionado tiene la misma estructura y el mismo orden de sumas que el
secuencial, y la igualdad vuelve a ser exacta. Cualquier interacción faltante, o
una partícula perdida en la migración, o una caja raíz que dejó de ser común,
aparece ahí como un error distinto de cero. Fue lo que dio confianza para
reescribir la migración en la semana 5.

Los tests 1-9 corren en el binario secuencial; los 10-19 requieren `mpirun`.

## Interfaz de línea de comandos

```
-n N                 partículas (default 1000)
-dt DT               paso temporal (0.001)
-t T_END             tiempo final (1.0)
-s EPS               softening (0.01)
-o FILE              salida CSV
--init TYPE          plummer | uniform | twobody
--method M           bh | direct   (direct solo secuencial: es el oráculo)
--theta T            ángulo de apertura (0.5)
--threads T          hilos OpenMP
--rebalance-every K  repartición cada K pasos (20)
--decomp D           orb | morton          [MPI]
--exchange E         let | replicate       [MPI]
--balance B          work | count          [MPI]
--theta-let T        θ del criterio de exportación (default: = --theta)
--metrics            desglose por fase     [MPI]
--validate           suite de tests
--seed S             semilla (42)
```

## Resultados por semana

| Semana | Qué se agregó | Hallazgo |
|---|---|---|
| 1 | Secuencial, leapfrog, O(N²), Morton | 4/4 tests |
| 2 | Barnes-Hut, octree | 8/8; speedup 14,8× vs directo a N=100K |
| 3 | MPI + OpenMP, réplica global | 13/13. **El cuello no era la comunicación** (Allgatherv 1,3%) sino el árbol replicado |
| 4 | LET, balance por costo | 17/18. **El LET era correcto pero no comprimía** (k=0,92): los dominios Morton no son compactos |
| 5 | ORB | 19/19. k de 0,92 a **0,66**; `tree_time` de crecer 4,8× a crecer 1,6× |

Dos veces seguidas la hipótesis razonable estaba equivocada, y las dos veces lo
mostró una medición barata: el desglose por fase en la semana 3, y en la semana 4
que `--theta-let` no cambiara nada más el volumen de las cajas de dominio.

El punto que vale la pena subrayar: **el LET de la semana 4 era correcto** —θ=0
bit a bit, masa a 5e-14, energía conservada— **y aun así no servía**. Correcto y
efectivo son propiedades distintas. La suite certificaba la primera y era ciega a
la segunda; lo que detectó el problema fue una métrica, no un test.
