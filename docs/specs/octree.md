# Módulo: octree

## Descripción general

Construye un octree (árbol de 8 hijos) que subdivide recursivamente el dominio cúbico hasta que cada hoja contiene a lo sumo una partícula. Cada nodo almacena la masa total y el centro de masa de su subárbol, información que Barnes-Hut usa para aproximar la fuerza de grupos lejanos por un único cuerpo equivalente. Es la estructura central de la etapa 4 del plan y la base espacial sobre la que opera el módulo `barnes_hut`.

## Conceptos

### Pool de nodos con índices

En lugar de nodos enlazados por punteros, se usa un arreglo dinámico `OctreeNode *nodes` e índices `int` (campos `children[8]`, `root`). Ventajas:

- Menos llamadas a `malloc`/`free` (una reasignación amortizada por duplicación).
- Mejor localidad de caché en el recorrido.
- Serialización directa para el intercambio LET por MPI (semanas 3-4).

**Precaución**: `realloc` puede reubicar `nodes`. El código nunca cachea punteros a `nodes[]` a través de una inserción; siempre reindexado por `t->nodes[idx]`.

### Estados de un nodo

| Estado        | `is_leaf` | `particle` | `children[]` |
|---------------|-----------|------------|--------------|
| Hoja vacía    | 1         | -1         | todos -1     |
| Hoja ocupada  | 1         | ≥0         | todos -1     |
| Nodo interno  | 0         | -1         | 8 índices    |

### Cadena de partículas (`next[]`)

En el caso normal cada hoja tiene una sola partícula. Para partículas **coincidentes** (o al alcanzar `MAX_DEPTH = 21`, el límite de resolución de Morton), varias partículas comparten una hoja formando una lista enlazada: `t->next[i]` es la siguiente partícula de la hoja cuya cabeza es `node->particle`, con `-1` como terminador. Esto evita la recursión infinita en configuraciones degeneradas.

### Convención de octantes

El hijo correspondiente a una posición se indexa por los bits de signo respecto al centro de la celda:

```
octante = (x ≥ cx ? 1 : 0) | (y ≥ cy ? 2 : 0) | (z ≥ cz ? 4 : 0)
```

## Interfaz pública

### `octree_bounds(p, n, min, max)`

Calcula un bounding box **cúbico** que contiene a las `n` partículas, centrado en el centro del AABB, con lado igual a la mayor extensión de los 3 ejes más un margen del 1% por lado (evita que una partícula caiga exactamente sobre `max[k]`). Un cubo simétrico garantiza subdivisión uniforme en octantes.

### `octree_build(p, n)`

Reserva un `Octree` e inserta las `n` partículas una a una en el cubo global (inserción incremental con subdivisión perezosa). No calcula masas; llamar a `octree_compute_mass` después. El caller libera con `octree_free`.

### `octree_compute_mass(t, p)`

Recorre el árbol en **post-orden**: cada hoja agrega la masa y el momento (masa·posición) de las partículas de su cadena; cada nodo interno combina la masa total y el centro de masa ponderado de sus 8 hijos. Debe llamarse tras `octree_build` y antes de calcular fuerzas.

### `octree_free(t)`

Libera el pool `nodes`, el arreglo `next` y el propio `Octree`. Seguro con `NULL`.

## Funciones internas

| Función | Descripción |
|---|---|
| `alloc_node(t)` | Reserva un nodo (creciendo el pool por duplicación) inicializado como hoja vacía. Retorna su índice. |
| `octant_of(t, idx, pos)` | Octante 0..7 de `pos` respecto al centro de la celda del nodo `idx`. |
| `subdivide(t, idx)` | Convierte una hoja en nodo interno creando sus 8 subcubos. Captura el bbox por valor porque `alloc_node` puede reubicar el pool. |
| `octree_insert(t, idx, p, i, depth)` | Inserta la partícula `i` descendiendo; subdivide al colisionar o encadena en `next[]` si `depth ≥ MAX_DEPTH`. |
| `compute_mass_rec(t, idx, p)` | Recursión post-orden de `octree_compute_mass`. |

## Estructuras de datos

```c
typedef struct {
    double min[3], max[3];   // bounding box cúbico
    double mass;             // masa del subárbol
    double cm[3];            // centro de masa del subárbol
    int    children[8];      // índices a hijos; -1 si no subdividido
    int    particle;         // cabeza de cadena o -1
    int    is_leaf;
} OctreeNode;

typedef struct {
    OctreeNode *nodes;
    int         count, capacity;
    int        *next;        // cadena de partículas por hoja (tam n)
    int         root;
    double      size;        // lado del cubo raíz
} Octree;
```

## Dependencias

- `particle.h` — tipo `Particle`.
- `vec3.h` — macros vectoriales (`VEC3_ZERO`, `VEC3_COPY`).
- `<stdlib.h>` — `malloc`, `realloc`, `free`.

## Consideraciones

- **Memoria**: la subdivisión crea siempre los 8 hijos aunque algunos queden vacíos; el consumo es ~O(N) nodos. Para N=100K son ~16 KB–decenas de MB según distribución. Optimizar (hijos perezosos) queda fuera de alcance de la semana 2.
- **Reconstrucción por paso**: en el loop de simulación el árbol se rearma completo cada paso porque las partículas se mueven. Es lo correcto y simple; medir el costo relativo construcción vs fuerza.
- **`MAX_DEPTH = 21`** coincide con la resolución de Morton (21 bits/eje); a esa profundidad las partículas restantes se encadenan en la misma hoja.
