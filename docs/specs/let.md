# Spec: LET (Locally Essential Tree)

Etapa 8. Reemplaza la réplica global de la semana 3 (`replicate_particles`).

## Idea

Cada proceso recibe solo el subconjunto del árbol global que necesita para
calcular las fuerzas de **sus** partículas con la misma precisión que tendría con
el árbol entero. Lo lejano llega resumido en pocos nodos agregados.

## Precondición: caja raíz compartida

`octree_build` deriva el cubo raíz de las partículas que recibe. Con replicación
todos los procesos pasaban el mismo arreglo global y obtenían el mismo cubo; con
LET cada uno construye sobre sus N/P locales y obtendría un cubo distinto, con lo
que las celdas no alinearían entre procesos.

Por eso existe `octree_build_box(p, n, min, max)`, que toma el cubo explícito.
Todos los procesos pasan `Domain::gmin/gmax`. `octree_build` queda como wrapper.

## Criterio de exportación (MAC conservador)

Para el nodo `c` de lado `s` y centro de masa `cm`, y la caja destino `B_r`:

```
d_min = distancia mínima de cm a B_r      (0 si cm cae dentro de B_r)
aceptar el resumen si   s² < theta² · (d_min² + eps²)
```

**Correcto porque** el receptor, al calcular la fuerza sobre su partícula
`i ∈ B_r`, aplica `s² < theta²·(d_i² + eps²)` con `d_i ≥ d_min`. Usar la distancia
mínima sobre toda la caja garantiza que nunca se exporta un resumen donde el
receptor quería detalle. El caso inverso (abrir un nodo que el receptor habría
aceptado) manda datos de más: cuesta ancho de banda, no precisión.

Con `theta = 0` la condición nunca se cumple, se desciende siempre hasta las hojas
y el LET degenera en réplica completa. Es lo que hace exacto al test 14.

## Recorrido

```
si nodo vacío (mass == 0):        podar
si hoja:                          exportar las partículas REALES
si interno y MAC(c, B_r):         exportar PSEUDO-PARTÍCULA (mass, cm)
si interno y no MAC:              descender a los 8 hijos
```

## Representación

Lo exportado viaja como `Particle`, así se reutilizan el tipo MPI derivado y
`MPI_Alltoallv` sin empaquetado a mano:

| Campo | Nodo resumido | Hoja |
|---|---|---|
| `pos` | centro de masa del nodo | posición real |
| `mass` | masa del subárbol | masa real |
| `id` | `PARTICLE_GHOST_ID` (-1) | id real |
| `vel`, `acc`, `work` | 0 | valores reales |

## Orden dentro del paso

```
kick, drift
[cada K pasos] repartición (bbox global + splitters)
migración
octree_build_box(local, n_local, gmin, gmax)
octree_compute_mass                     <- el MAC de exportación lo necesita
let_exchange                            <- fantasmas appendeados en [n_local, n_local+n_ghost)
octree_insert_particles                 <- sin reconstruir el árbol
octree_compute_mass                     <- las masas cambiaron
compute_forces_bh_range(..., 0, n_local)
kick
```

Los fantasmas se insertan en el árbol ya construido en vez de armar uno fusionado:
la caja raíz es la misma, así que la estructura sigue válida y solo se subdividen
las hojas donde caen los puntos nuevos.

## Los fantasmas son transitorios

El `migrate_particles` del paso siguiente reemplaza el arreglo local y los
descarta sin limpieza explícita. Todo lo que sea **checksum, energía o migración
debe usar el rango `[0, n_local)`**: un fantasma tiene `id = -1` y la masa de un
nodo agregado, así que contarlo rompe el checksum e infla la energía.

## Se pierde la exactitud bit a bit

Un fantasma llega como punto de masa en `cm`; al reinsertarlo puede quedar
agrupado con otros bajo un ancestro que a su vez satisface el MAC, y aproximarse
dos veces. El error sigue acotado por θ pero ya no es cero, así que el criterio
del test 13 (error exactamente 0.0) no se puede sostener. Escalera de validación:

| Nivel | Compara | Criterio |
|---|---|---|
| θ=0 (test 14) | LET vs BH secuencial | **bit a bit** |
| θ=0.5 (test 15) | LET vs BH secuencial, y ambos vs directo | < 1e-3, y < 1.5× el error del BH |
| dinámico (test 18) | energía tras 200 pasos | drift < 1% |

## Verificación barata

`let_mass_error`: Σ mass sobre locales + fantasmas debe dar la masa total del
sistema, porque el LET resume pero no descarta. Detecta casi cualquier bug de
selección sin comparar contra una referencia secuencial. Medido: **4.9e-14**.

## Historia: el LET dependía de la descomposición

Con rangos Morton (semana 4) este mismo código **no comprimía**: `n_ghost` crecía
como N^0,92, la misma asintótica que la réplica, porque un tramo Morton contiguo
no es una región compacta y el AABB del destino cubría el 44,6% del dominio.

La semana 5 reemplazó la descomposición por [ORB](orb.md) y `let.c` **no necesitó
un solo cambio**: el criterio ya era correcto, lo que fallaba era la geometría que
recibía. Con dominios que son cajas, el exponente bajó a **0,66**, cerca del 2/3
de un efecto de superficie.

Es el ejemplo más claro del trabajo de que correcto y efectivo son propiedades
distintas: la suite certificaba lo primero (θ=0 bit a bit, masa a 5e-14) y era
completamente ciega a lo segundo.
