# Módulo: morton

## Descripción general

Implementa la codificación Morton (Z-order curve) para mapear coordenadas 3D a claves enteras de 63 bits que preservan la localidad espacial. El ordenamiento de partículas por clave Morton mejora la coherencia de caché durante el cálculo de fuerzas y es la base para la construcción eficiente de octrees (necesarios en Barnes-Hut).

## Conceptos

### Curva Z-order (Morton)

La curva de Morton recorre el espacio 3D intercalando los bits de las coordenadas discretizadas (x, y, z). Para 21 bits por eje, el patrón de bits de la clave es:

```
bit 62 ... bit 2  bit 1  bit 0
x₂₀ y₂₀ z₂₀ ... x₀    y₀    z₀
```

Partículas cercanas en el espacio 3D tienden a tener claves Morton cercanas, lo que al ordenar el arreglo por estas claves agrupa partículas vecinas en posiciones contiguas de memoria.

### Bit interleaving

La función interna `expand_bits` dispersa los 21 bits inferiores de un entero en las posiciones 0, 3, 6, 9, ..., 60 usando una secuencia de operaciones de shift y mask. Luego, la clave final se construye como:

```
morton = (expand(ix) << 2) | (expand(iy) << 1) | expand(iz)
```

## Interfaz pública

### `morton_encode(x, y, z, min, max, bits)`

Codifica una posición 3D continua en una clave Morton entera.

**Proceso**:
1. Normaliza `(x, y, z)` al rango `[0, 2^bits - 1]` usando el bounding box `[min, max]`.
2. Convierte a enteros `(ix, iy, iz)`.
3. Intercala los bits de los 3 ejes para producir una clave de `3 * bits` bits.

**Parámetros**:
- `x, y, z`: coordenadas del punto en espacio continuo
- `min[3], max[3]`: esquinas del bounding box que define el dominio
- `bits`: resolución por eje (típicamente 21, produciendo claves de 63 bits)

**Retorna**: clave Morton `uint64_t`.

### `morton_sort(p, n, min, max)`

Asigna claves Morton (con 21 bits por eje) a todas las `n` partículas según su posición y luego ordena el arreglo in-place por clave creciente usando `qsort`.

**Parámetros**:
- `p`: arreglo de partículas (lectura de `pos[]`, escritura de `morton`, reordenamiento completo)
- `n`: cantidad de partículas
- `min[3], max[3]`: bounding box del dominio

## Funciones internas

| Función | Descripción |
|---|---|
| `expand_bits(v)` | Dispersa los 21 bits inferiores de `v` en posiciones 0, 3, 6, ..., 60 mediante 5 etapas de shift+mask. Complejidad O(1). |
| `morton_compare(a, b)` | Comparador para `qsort` que ordena `Particle` por campo `morton` ascendente. |

## Dependencias

- `particle.h` — tipo `Particle`.
- `<stdint.h>` — `uint64_t`.
- `<stdlib.h>` — `qsort`.

## Consideraciones

- La resolución de 21 bits por eje da 2²¹ = 2.097.152 celdas por dimensión, suficiente para discriminar partículas en simulaciones de hasta ~10⁷ cuerpos.
- El bounding box debe recalcularse cada vez que las partículas se mueven significativamente; en el loop principal se calcula una vez al inicio.
- Partículas exactamente en `max[k]` se mapean al bin máximo gracias a que el rango se calcula como `2^bits - 1`.
