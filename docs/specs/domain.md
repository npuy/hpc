# Módulo: domain

## Descripción general

Reparte el espacio de claves Morton `[0, 2^63)` entre los procesos MPI, un tramo contiguo por proceso. Es la etapa 6 del plan general (MPI estático) y la base sobre la que operan `migration` (etapa 7) y el rebalanceo dinámico de la semana 4 (etapa 9).

## Conceptos

### Por qué tramos contiguos de clave Morton

La curva Z preserva localidad: puntos cercanos en el espacio 3D tienden a tener claves cercanas. Un tramo contiguo de claves es entonces una región compacta del espacio, no un conjunto disperso de celdas. Esto es lo que hace baratos la migración (una partícula que se mueve poco cambia poco su clave, y casi nunca cruza de tramo) y el LET de la semana 4 (los vecinos de un proceso son pocos procesos, no todos).

### Reparto por conteo, no por volumen

Cortar el espacio de claves en `P` partes iguales sería trivial y estaría mal. Una distribución de Plummer concentra casi todas las partículas en unas pocas celdas centrales: el proceso que recibiera el centro tendría casi todo el trabajo y el resto quedaría ocioso.

`domain_partition` busca en cambio los `P-1` splitters que dejan `N/P` partículas en cada tramo. Medido en el proyecto: con Plummer y `P ∈ {4, 8}` el desbalance `max/avg` resultante es **1.0000**.

### Bisección por histograma

El algoritmo busca, para cada `k`, la menor clave `S` tal que la cantidad de partículas con clave `< S` alcanza el objetivo `k·N/P`. Es una búsqueda binaria sobre el espacio de claves donde el "test" es un conteo global:

```
lo = 0, hi = 2^63
mientras lo < hi:
    mid = (lo + hi) / 2
    c_local  = # particulas locales con clave < mid     (busqueda binaria, O(log n))
    c_global = MPI_Allreduce(c_local, SUM)
    si c_global >= objetivo: hi = mid
    si no:                   lo = mid + 1
```

Converge en 63 iteraciones (una por bit de la clave).

**Los `P-1` splitters se bisectan simultáneamente**: cada iteración hace un solo `MPI_Allreduce` de `P-1` enteros en vez de `P-1` reducciones separadas. El costo total es de ~63 colectivos, independiente de `P`.

Ventajas frente a un sample sort: no requiere ordenar globalmente ni juntar las claves en un proceso, y es exactamente el procedimiento que reusará el rebalanceo dinámico — en la semana 4 basta con volver a llamarlo cada K pasos.

### Terminación colectiva sin divergencia

`lo[]` y `hi[]` se derivan únicamente de valores globales (`c_global`, idéntico en todos los procesos por ser resultado de un `Allreduce`). Por lo tanto todos los procesos recorren exactamente las mismas iteraciones y evalúan la misma condición de corte. Esto es lo que garantiza que ningún proceso abandone el bucle colectivo antes que otro, que sería un deadlock.

### Bounding box congelado

`domain_global_bounds` calcula el AABB global con `MPI_Allreduce(MIN/MAX)` y lo convierte en cubo con margen del 1%, con la misma construcción que `octree_bounds` para que árbol y particionamiento vean la misma geometría.

El box se **congela** entre reparticiones. Las claves solo son comparables entre procesos si todos usan el mismo box; si se recodificaran las claves con un box distinto al que produjo los splitters, la asignación de dueños sería inconsistente y se perderían partículas. Por eso `repartition()` en `main.c` actualiza siempre box y splitters **juntos**.

La contrapartida es que entre reparticiones una partícula puede salirse del box congelado. `morton_encode` acota las coordenadas al rango válido antes de discretizar; sin ese clamp, una partícula fuera del box recibiría una clave envuelta (extremo derecho → clave del extremo izquierdo) y migraría al proceso equivocado.

## API

| Función | Rol |
|---|---|
| `domain_init` | Reserva `splitters[nprocs+1]`, consulta rank/nprocs |
| `domain_free` | Libera `splitters` |
| `domain_global_bounds` | AABB global cúbico por `Allreduce`; congela `gmin`/`gmax` |
| `domain_partition` | Bisección por histograma; llena `splitters` |
| `domain_owner` | Dueño de una clave; búsqueda binaria O(log P) |

## Estructura `Domain`

| Campo | Significado |
|---|---|
| `rank`, `nprocs` | Identidad dentro del comunicador |
| `gmin[3]`, `gmax[3]` | Bounding box global cúbico congelado |
| `splitters` | `nprocs+1` cortes; el tramo de `k` es `[splitters[k], splitters[k+1])` |

`splitters[0] == 0` y `splitters[nprocs] == UINT64_MAX` son centinelas: garantizan que toda clave tenga dueño.

## Precondiciones

- `domain_partition` asume el arreglo local **ordenado por clave creciente** (usa búsqueda binaria para contar) y las claves calculadas con `gmin`/`gmax`.
- Ambas funciones son **colectivas**: todos los procesos del comunicador deben invocarlas.

## Casos borde

| Caso | Comportamiento |
|---|---|
| `nprocs == 1` | `domain_partition` retorna de inmediato; el único tramo es todo el espacio |
| `n_local == 0` | Legal. Aporta los neutros de MIN/MAX y cuenta 0 en las bisecciones |
| Claves repetidas | Dos splitters pueden coincidir ⇒ tramo vacío. Los splitters son **no decrecientes**, no estrictamente crecientes |
| Más procesos que partículas | Verificado con `N=3, P=4`: los procesos sobrantes quedan con 0 partículas |
