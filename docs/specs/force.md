# Módulo: force

## Descripción general

Implementa el cálculo de interacciones gravitatorias entre partículas usando el método directo (all-pairs) con complejidad O(N²). Aprovecha la tercera ley de Newton (simetría acción-reacción) para reducir el número de evaluaciones a N(N-1)/2. También provee funciones para calcular las energías cinética y potencial del sistema, necesarias para la validación de conservación de energía.

## Interfaz pública

### `compute_forces_direct(p, n, softening)`

Calcula la aceleración gravitatoria sobre cada partícula por interacción directa con todas las demás. Usa softening gravitacional para evitar singularidades cuando dos partículas están muy cerca.

**Algoritmo**:
1. Resetea `acc[3]` de todas las partículas a cero.
2. Para cada par `(i, j)` con `j > i`:
   - Calcula el vector desplazamiento `r = pos_j - pos_i`
   - Calcula la distancia softened: `dist² = |r|² + ε²`
   - Calcula `1/dist³` para obtener la magnitud de la fuerza
   - Acumula `acc_i += m_j / dist³ * r` y `acc_j -= m_i / dist³ * r`

**Softening**: El parámetro `ε` (softening) se eleva al cuadrado y se suma a `|r|²` antes de calcular la inversa. Esto regulariza la fuerza gravitatoria en distancias cortas, previniendo aceleraciones divergentes y mejorando la estabilidad numérica del integrador.

**Complejidad**: O(N²) en tiempo, O(1) en espacio auxiliar.

### `kinetic_energy(p, n)`

Calcula la energía cinética total del sistema: `K = Σ ½ mᵢ |vᵢ|²`.

### `potential_energy(p, n, softening)`

Calcula la energía potencial gravitatoria total: `U = -Σᵢ<ⱼ mᵢ mⱼ / √(|rᵢⱼ|² + ε²)`.

Usa el mismo softening que `compute_forces_direct` para que la energía total `E = K + U` sea una integral de movimiento consistente con las ecuaciones integradas.

## Dependencias

- `particle.h` — tipo `Particle`.
- `vec3.h` — macros `VEC3_SUB`, `VEC3_DOT`, `VEC3_ZERO`.
- `<math.h>` — `sqrt`.

## Consideraciones de rendimiento

- El doble bucle con simetría de Newton reduce las evaluaciones de fuerza a la mitad respecto a un bucle naïve.
- Es el bottleneck computacional de la simulación: para N=1000 y 1000 pasos, se evalúan ~500 millones de interacciones.
- En la versión paralela, este módulo será reemplazado/complementado por el algoritmo Barnes-Hut O(N log N).
