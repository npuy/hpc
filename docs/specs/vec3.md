# Módulo: vec3

## Descripción general

Librería header-only de aritmética vectorial 3D implementada con macros del preprocesador. Proporciona las operaciones básicas necesarias para la simulación N-body: suma, resta, escalado, multiply-add, producto punto y norma euclidiana.

## Justificación de diseño

Se usan macros en lugar de funciones inline para garantizar zero-overhead en compiladores sin soporte completo de `inline`, y para que las operaciones puedan aplicarse sobre cualquier arreglo `double[3]` sin necesidad de un tipo opaco. El patrón `do { ... } while(0)` asegura que las macros se comporten como sentencias únicas en contextos de control de flujo (`if`, `for`, etc.).

## Interfaz pública

| Macro | Firma lógica | Descripción |
|---|---|---|
| `VEC3_ZERO(a)` | `a → void` | Pone las 3 componentes de `a` en cero |
| `VEC3_COPY(r, a)` | `(r, a) → void` | Copia `a` en `r` componente a componente |
| `VEC3_ADD(r, a, b)` | `(r, a, b) → void` | `r = a + b` |
| `VEC3_SUB(r, a, b)` | `(r, a, b) → void` | `r = a - b` |
| `VEC3_SCALE(r, a, s)` | `(r, a, s) → void` | `r = a * s` donde `s` es escalar |
| `VEC3_MADD(r, a, b, s)` | `(r, a, b, s) → void` | `r = a + b * s` (multiply-add) |
| `VEC3_DOT(a, b)` | `(a, b) → double` | Producto punto `a · b` |
| `VEC3_NORM(a)` | `a → double` | Norma euclidiana `‖a‖ = sqrt(a · a)` |

## Dependencias

- `<math.h>` — para `sqrt()` en `VEC3_NORM`.

## Notas de uso

- Los argumentos deben ser arreglos `double[3]` o punteros a `double` con al menos 3 elementos.
- Las macros evalúan sus argumentos múltiples veces; no pasar expresiones con efectos secundarios (ej: `VEC3_ADD(r, a++, b)` tiene comportamiento indefinido).
- `VEC3_DOT` y `VEC3_NORM` son expresiones (no sentencias), por lo que pueden usarse como parte de asignaciones o condicionales.
