# Semana 2 — Reporte: Barnes-Hut Secuencial

## Archivos implementados

```
src/
  octree.h / .c       -- octree (pool de nodos, bounding box cubico, centros de masa)
  barnes_hut.h / .c   -- fuerza O(N log N) por traversal con criterio de apertura (MAC)
  validation.h / .c   -- +4 tests (5-8): masa, centro de masa, error BH, convergencia theta
  main.c              -- flags --method direct|bh y --theta; arbol reconstruido por paso
Makefile              -- octree.c y barnes_hut.c agregados a la compilacion
docs/specs/
  octree.md           -- spec del modulo octree
  barnes_hut.md       -- spec del modulo barnes_hut
```

Se conserva `compute_forces_direct` (O(N²)) intacto como **oráculo** de referencia.

## Resultados de validacion

| Test | Descripcion | Criterio | Resultado | Valor obtenido |
|------|-------------|----------|-----------|----------------|
| 1 | Orbita circular 2 cuerpos | Error posicion < 1% | **PASS** | 0.001% |
| 2 | Conservacion de energia (N=100) | Drift < 0.1% | **PASS** | 3.12e-07 |
| 3 | Conservacion de momento lineal | Drift < 1e-12 | **PASS** | 3.54e-17 |
| 4 | Morton ordering (8 esquinas) | Orden correcto | **PASS** | 8/8 Z-order |
| 5 | Masa del arbol == masa total | Error < 1e-12 | **PASS** | 4.44e-16 |
| 6 | Centro de masa del arbol | ‖·‖ < 1e-10 | **PASS** | 1.06e-16 |
| 7 | Error fuerza BH vs directo (θ=0.5) | Error L2 < 1% | **PASS** | 3.98e-03 (0.40%) |
| 8 | Convergencia en θ | Monotono + θ=0 exacto | **PASS** | ver abajo |

**8/8 tests pasaron.**

### Test 8: convergencia en θ

| θ    | Error L2 vs directo |
|------|---------------------|
| 0.8  | 1.64e-02            |
| 0.5  | 3.98e-03            |
| 0.2  | 2.01e-04            |
| 0.0  | 1.15e-15            |

El error decrece monótonamente al bajar θ y **θ=0 reproduce el método directo a precisión de máquina** (1.15e-15), confirmando que la aproximación es consistente.

## Rendimiento: directo O(N²) vs Barnes-Hut O(N log N)

Tiempo por paso (ms), θ=0.5, distribución Plummer:

| N       | Directo (ms/paso) | Barnes-Hut (ms/paso) | Speedup |
|---------|-------------------|----------------------|---------|
| 1,000   | 2.03              | 4.24                 | 0.47×   |
| 5,000   | 48.8              | 38.7                 | 1.26×   |
| 10,000  | 197.4             | 95.4                 | 2.07×   |
| 50,000  | 5,057             | 640                  | 7.90×   |
| 100,000 | 22,744            | 1,538                | 14.8×   |

- **Punto de cruce**: entre N=1K y N=5K. Para N pequeño el costo de construir el árbol domina y Barnes-Hut es más lento; a partir de ~N=5K ya conviene.
- **Escalamiento directo**: 50K→100K multiplica el tiempo por ~4.5× (esperado ~4× por O(N²)).
- **Escalamiento BH**: 10K→100K (10× partículas) multiplica el tiempo por ~16×, consistente con O(N log N) (muy por debajo del 100× de O(N²)).

## Conservacion de energia con Barnes-Hut

| Config | Pasos | Energy drift |
|--------|-------|--------------|
| N=2,000, θ=0.5, dt=0.001 | 500 | 7.88e-05 |
| N=100,000, θ=0.5, dt=0.001 | 5 | 1.12e-07 |

La aproximación de Barnes-Hut introduce ruido acotado; el integrador leapfrog mantiene el drift pequeño (< 0.01%). El árbol se reconstruye cada paso.

## Verificacion de memoria

`leaks` (macOS) sobre un run BH completo: **0 leaks for 0 total leaked bytes**. El árbol se libera correctamente con `octree_free` en cada paso.

## Componentes implementados

### Octree (`octree.c`)
- Pool de nodos con índices `int` (no punteros): serializable para LET/MPI, robusto ante `realloc`.
- Bounding box cúbico centrado con margen del 1%.
- Inserción incremental con subdivisión perezosa; cadena `next[]` para partículas coincidentes o profundidad máxima (`MAX_DEPTH=21`).
- Centros de masa por recorrido post-orden.

### Barnes-Hut (`barnes_hut.c`)
- Traversal con criterio de apertura `s² < θ²·d²` (sin raíces).
- Misma fórmula con softening que el método directo ⇒ θ=0 lo reproduce exactamente.
- `potential_energy_bh` para monitoreo de energía en O(N log N).

### Integracion (`main.c`)
- Flags `--method direct|bh` (default `bh`) y `--theta` (default 0.5).
- Árbol construido/liberado cada paso; energía calculada con el método activo.

## Criterios de aceptacion — Resumen

| Criterio | Estado |
|----------|--------|
| Masa de la raíz == masa total (< 1e-12) | OK (4.44e-16) |
| CM de la raíz == CM directo (< 1e-10) | OK (1.06e-16) |
| Error fuerza BH vs directo, θ=0.5 (< 1%) | OK (0.40%) |
| θ→0 reproduce el método directo (< 1e-9) | OK (1.15e-15) |
| Error decrece monótonamente al bajar θ | OK |
| Conservación de energía con BH (< 1%) | OK (7.88e-05) |
| BH más rápido que directo para N≥50K | OK (7.9×–14.8×) |
| N=100K con BH completa sin crash | OK |
| Sin fugas de memoria | OK (leaks limpio) |

**Todos los criterios de aceptacion se cumplen.**

## Argumentos CLI (nuevos)

```
--method M    Metodo de fuerza: direct|bh (default: bh)
--theta T     Angulo de apertura Barnes-Hut (default: 0.5)
```

## Siguiente paso: Semana 3

Con Barnes-Hut secuencial validado, la semana 3 aborda la **versión híbrida básica**:
paralelización intra-nodo con OpenMP (recorrido del árbol por partícula es
trivially parallel) y distribución de partículas con MPI (particionamiento por rangos
de clave Morton + migración). El pool de nodos con índices ya está preparado para la
serialización que requerirá el intercambio LET.
