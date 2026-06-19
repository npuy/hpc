# Módulo: validation

## Descripción general

Suite de tests de validación física y numérica para verificar la correctitud de la simulación N-body. Cada test verifica una propiedad conservada o un comportamiento esperado del sistema, sirviendo como regression tests ante cambios en los módulos de fuerza, integración o inicialización.

## Tests

### Test 1: Órbita circular de dos cuerpos (`test_two_body_orbit`)

**Qué verifica**: Que dos cuerpos de igual masa, inicializados en órbita circular, regresen a su posición inicial después de un período orbital completo.

**Metodología**:
1. Inicializa el sistema con `init_two_body` (2 masas iguales, separación unitaria).
2. Calcula el período orbital analítico: `T = 2π r^(3/2) / √(GM)`.
3. Integra 10.000 pasos con Leapfrog KDK.
4. Mide el error de posición de ambos cuerpos respecto a la posición inicial.

**Criterio de aceptación**: Error relativo de posición < 1% de la separación inicial.

**Softening usado**: ε = 0.001 (pequeño para no distorsionar la órbita).

### Test 2: Conservación de energía (`test_energy_conservation`)

**Qué verifica**: Que la energía total `E = K + U` no presente drift secular durante la integración.

**Metodología**:
1. Genera 100 partículas con distribución de Plummer (seed=42).
2. Integra 100 pasos con dt=0.001, softening=0.01.
3. Compara la energía final `Ef` con la inicial `E0`.

**Criterio de aceptación**: `|Ef - E0| / |E0| < 10⁻³`.

**Fundamento**: El integrador Leapfrog es simpléctico, por lo que la energía debe oscilar sin drift secular. El umbral de 10⁻³ es conservador; con dt más pequeño el drift sería menor.

### Test 3: Conservación de momento lineal (`test_momentum_conservation`)

**Qué verifica**: Que el momento lineal total `P = Σ mᵢvᵢ` se conserve exactamente (a precisión de máquina).

**Metodología**:
1. Genera 50 partículas con Plummer (seed=123).
2. Integra 50 pasos con dt=0.001, softening=0.01.
3. Compara el momento total final con el inicial.

**Criterio de aceptación**: `|Pf - P0| < 10⁻¹²`.

**Fundamento**: La simetría de Newton en `compute_forces_direct` (acción = -reacción) garantiza conservación exacta del momento. Cualquier error por encima de la precisión de máquina indica un bug en el cálculo de fuerzas.

### Test 4: Ordenamiento Morton (`test_morton_ordering`)

**Qué verifica**: Que la codificación Morton asigne claves distintas y correctamente ordenadas a las 8 esquinas del cubo unitario.

**Metodología**:
1. Crea 8 partículas en las esquinas de `[0,1]³`.
2. Calcula claves Morton con 21 bits por eje.
3. Aplica `morton_sort`.
4. Verifica que las claves estén en orden estrictamente creciente y que las 8 sean distintas.

**Criterio de aceptación**: Orden estricto y 8 claves únicas.

### Runner: `run_all_tests`

Ejecuta los 4 tests secuencialmente e imprime un resumen `pasados/total`. Retorna 0 si todos pasan, 1 si alguno falla. Se invoca desde `main` con la flag `--validate`.

## Dependencias

Usa todos los módulos del proyecto: `particle`, `force`, `leapfrog`, `morton`, `vec3`.

## Cómo ejecutar

```bash
make
./nbody --validate
```
