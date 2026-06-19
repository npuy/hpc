# Módulo: main

## Descripción general

Punto de entrada de la simulación secuencial de N cuerpos. Parsea los argumentos de línea de comandos, inicializa las partículas según la distribución elegida, ejecuta el loop de integración temporal con el esquema Leapfrog KDK, y reporta métricas de rendimiento y conservación de energía.

## Argumentos de línea de comandos

| Flag | Tipo | Default | Descripción |
|---|---|---|---|
| `-n N` | int | 1000 | Número de partículas |
| `-dt DT` | double | 0.001 | Paso temporal |
| `-t T_END` | double | 1.0 | Tiempo final de simulación |
| `-s EPS` | double | 0.01 | Parámetro de softening gravitacional |
| `-o FILE` | string | (ninguno) | Prefijo para archivos de salida CSV |
| `--init TYPE` | string | plummer | Distribución inicial: `plummer`, `uniform`, `twobody` |
| `--validate` | flag | off | Ejecuta la suite de validación en vez de simular |
| `--seed SEED` | int | 42 | Semilla para el generador de números aleatorios |
| `-h, --help` | flag | — | Muestra el mensaje de uso |

## Flujo de ejecución

1. **Parseo de argumentos**: Lee flags de la CLI con un loop manual sobre `argv`.
2. **Modo validación**: Si `--validate` está presente, delega a `run_all_tests()` y termina.
3. **Inicialización de partículas**: Según `--init`, llama a `init_plummer`, `init_uniform_cube` o `init_two_body`.
4. **Ordenamiento Morton**: Calcula el bounding box del sistema y ordena las partículas por clave Morton para mejorar localidad de caché.
5. **Cálculo de fuerzas inicial**: `compute_forces_direct` para obtener las aceleraciones en t=0.
6. **Snapshot inicial** (si `-o`): Escribe el estado inicial a CSV.
7. **Loop de integración**:
   - Para cada paso: Kick(½dt) → Drift(dt) → Forces → Kick(½dt)
   - Cada `snapshot_interval` pasos (=100): imprime energía y drift relativo.
8. **Medición de tiempo**: Usa `clock_gettime(CLOCK_MONOTONIC)` para wall-clock time.
9. **Reporte final**: Imprime energía final, drift, tiempo total y tiempo por paso.
10. **Snapshot final** (si `-o`): Escribe el estado final a `{FILE}_final.csv`.

## Salida por consola

```
N-Body Simulacion Secuencial
  N = 1000, dt = 0.001, t_end = 1.0, softening = 0.01
  Init: plummer, seed: 42

Energia inicial E0 = -2.4567890123e-01
  paso 0/1000  t=0.0010  E=-2.4567890100e-01  drift=9.3456e-09
  paso 100/1000  t=0.1010  E=-2.4567891234e-01  drift=4.5234e-07
  ...

=== Resultado final ===
  Pasos: 1000
  Energia final Ef = -2.4567892345e-01
  Energy drift: 9.0567e-07
  Tiempo de ejecucion: 12.345 s
  Tiempo por paso: 12.345 ms
```

## Dependencias

Todos los módulos: `particle`, `force`, `leapfrog`, `morton`, `validation`, `vec3`.
