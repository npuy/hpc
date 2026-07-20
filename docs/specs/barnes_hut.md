# Módulo: barnes_hut

## Descripción general

Calcula las fuerzas gravitatorias en O(N log N) recorriendo el octree del módulo `octree`. Para cada partícula desciende desde la raíz aplicando el criterio de apertura: las celdas suficientemente lejanas se aproximan por su masa y centro de masa (un único cuerpo), y solo las cercanas se expanden a sus hijos. Implementa la etapa 2 del plan y reemplaza a `compute_forces_direct` (O(N²)) en el loop de simulación.

## Conceptos

### Criterio de apertura (MAC)

Para un nodo interno de lado `s` cuyo centro de masa está a distancia `d` de la partícula, se compara:

```
s / d < θ        (Multipole Acceptance Criterion)
```

- Si se cumple, la celda es **lo bastante lejana**: se acumula la fuerza usando `(mass, cm)` del nodo como un solo cuerpo.
- Si no, se **desciende** a los 8 hijos.

Para evitar la raíz cuadrada, se compara la forma equivalente `s² < θ²·d²`. La distancia usa `d² = |cm − pos|² + ε²` (mismo softening que la fuerza).

### Parámetro θ (theta)

| θ      | Efecto                                              |
|--------|-----------------------------------------------------|
| 0      | No aproxima nada ⇒ idéntico al método directo O(N²) |
| 0.5    | Estándar: buen balance error/velocidad              |
| > 0.5  | Más rápido, menos preciso                            |

Con θ=0 el recorrido llega siempre a las hojas y reproduce el directo salvo orden de sumas (verificado: error L2 ~10⁻¹⁵).

### Consistencia con el método directo

Se usa **exactamente la misma fórmula** con softening que `force.c`:

```
a_i += m_nodo · r / (|r|² + ε²)^{3/2},   r = cm_nodo − pos_i
```

Esto garantiza que el límite θ→0 coincida numéricamente con `compute_forces_direct` y que la comparación de error sea significativa.

## Interfaz pública

### `compute_forces_bh(t, p, n, theta, softening)`

Calcula la aceleración de cada partícula recorriendo el octree `t`. Resetea `acc[]` antes de acumular. El árbol debe tener masas/centros de masa ya calculados (`octree_compute_mass`).

**Parámetros**:
- `t`: octree construido (const, no se modifica).
- `p`, `n`: arreglo de partículas (escribe `acc[]`).
- `theta`: ángulo de apertura.
- `softening`: ε de regularización.

### `potential_energy_bh(t, p, n, theta, softening)`

Energía potencial total aproximada por el mismo recorrido: `U = ½ Σᵢ mᵢ φᵢ`, con `φᵢ` el potencial en la partícula `i`. Permite monitorear conservación de energía en O(N log N) en vez del O(N²) de `potential_energy`. Con θ=0 coincide con el directo.

## Funciones internas

| Función | Descripción |
|---|---|
| `accumulate_force(t, idx, i, p, θ², ε², acc)` | Recursión del cálculo de aceleración: hoja ⇒ suma directa sobre la cadena (salteando `i`); interno ⇒ MAC. |
| `potential_at(t, idx, i, p, θ², ε²)` | Recursión análoga que retorna el potencial en la partícula `i`. |

El árbol es inmutable durante el cálculo, por lo que cachear el puntero al nodo (`const OctreeNode *nd`) es seguro (no hay `realloc`).

## Dependencias

- `particle.h` — tipo `Particle`.
- `octree.h` — tipo `Octree`/`OctreeNode` y cadena `next[]`.
- `vec3.h` — macros vectoriales.
- `<math.h>` — `sqrt`.

## Complejidad

- **Tiempo**: O(N log N) promedio (cada partícula visita ~log N celdas con θ>0).
- **Espacio**: O(1) adicional (recursión de profundidad ≤ `MAX_DEPTH`).

## Consideraciones

- **Reconstrucción por paso**: el árbol se rearma cada paso en el loop principal; `compute_forces_bh` solo lo recorre.
- **Error controlado por θ**: verificado monótono decreciente al bajar θ (tests 7 y 8 de la suite).
- **Softening**: idéntico a `force.c` para consistencia del hamiltoniano; `potential_energy_bh` usa el mismo ε.
