# Semana 1 — Reporte: Version Secuencial

## Archivos implementados

```
src/
  vec3.h              -- macros inline para operaciones vectoriales 3D
  particle.h / .c     -- struct Particle, inicializacion (plummer, uniform, twobody), I/O CSV
  force.h / .c        -- fuerza directa O(N^2) con simetria Newton, energia cinetica/potencial
  leapfrog.h / .c     -- integrador leapfrog kick-drift-kick (KDK)
  morton.h / .c        -- codificacion Morton 3D (21 bits/eje), sort por Z-order
  validation.h / .c   -- suite de 4 tests de validacion
  main.c              -- programa principal con CLI y loop de simulacion
Makefile              -- targets: nbody (gcc), nbody_mpi (mpicc), clean
```

## Resultados de validacion

| Test | Descripcion | Criterio | Resultado | Valor obtenido |
|------|-------------|----------|-----------|----------------|
| 1 | Orbita circular 2 cuerpos (1 periodo) | Error posicion < 1% | **PASS** | 0.001% |
| 2 | Conservacion de energia (N=100, 100 pasos) | Drift < 0.1% | **PASS** | 3.12e-07 (0.00003%) |
| 3 | Conservacion de momento lineal | Drift < 1e-12 | **PASS** | 3.54e-17 |
| 4 | Morton ordering (8 esquinas cubo) | Orden correcto, claves unicas | **PASS** | 8/8 claves Z-order correctas |

**4/4 tests pasaron.**

## Rendimiento

| N | dt | t_end | Pasos | Tiempo total | Tiempo/paso | Energy drift |
|---|------|-------|-------|-------------|-------------|--------------|
| 1,000 | 0.001 | 0.1 | 100 | 0.194 s | 1.94 ms | 3.80e-08 |
| 1,000 | 0.001 | 1.0 | 1000 | 1.752 s | 1.75 ms | 1.40e-08 |
| 10,000 | 0.001 | 0.01 | 10 | 1.938 s | 193.8 ms | 2.60e-08 |

- **Escalamiento O(N^2):** N=1000 -> 1.75 ms/paso; N=10000 -> 193.8 ms/paso. Ratio ~110x para 10x particulas (esperado: ~100x). Consistente con O(N^2).
- **N=10000 < 30s/paso:** cumple ampliamente el criterio (193.8 ms << 30 s).

## Componentes implementados

### Particulas (`particle.c`)
- **init_two_body:** dos masas iguales en orbita circular, separacion unitaria.
- **init_plummer:** modelo de Plummer con rejection sampling para velocidades, correccion de centro de masa.
- **init_uniform_cube:** distribucion uniforme en cubo [-0.5, 0.5]^3.
- **I/O CSV:** escritura y lectura de snapshots con formato `t,mass,x,y,z,vx,vy,vz`.

### Fuerza directa (`force.c`)
- Calculo O(N^2) con simetria de Newton (N(N-1)/2 interacciones).
- Softening gravitacional para evitar singularidades.
- G = 1 (unidades naturales).
- Funciones de energia cinetica y potencial para validacion.

### Leapfrog KDK (`leapfrog.c`)
- Integrador simplectico kick-drift-kick.
- Conserva energia a largo plazo (drift ~1e-8 en 1000 pasos).
- Posiciones y velocidades sincronizadas al final de cada paso.

### Morton 3D (`morton.c`)
- Codificacion con 21 bits por eje (63 bits totales en uint64).
- Bit-interleaving con magic numbers para expand_bits.
- Sort de particulas por clave Morton (qsort).

### Validacion (`validation.c`)
- 4 tests automatizados cubriendo correctitud fisica y ordenamiento espacial.
- Ejecutable via `./nbody --validate`.

## Criterios de aceptacion — Resumen

| Criterio | Estado |
|----------|--------|
| Dos cuerpos completan orbita (error < 1%) | OK (0.001%) |
| Conservacion energia N=1000, 1000 pasos (drift < 0.1%) | OK (1.40e-08) |
| Conservacion momento lineal (drift < 1e-12) | OK (3.54e-17) |
| Morton: 8 esquinas del cubo en orden correcto | OK |
| Morton sort + verificacion | OK |
| Simulacion N=10000 completa sin crash | OK |
| Tiempo razonable N=10000 (< 30s/paso) | OK (0.194 s) |

**Todos los criterios de aceptacion se cumplen.**

## Argumentos CLI

```
./nbody [opciones]
  -n N          Numero de particulas (default: 1000)
  -dt DT        Paso temporal (default: 0.001)
  -t T_END      Tiempo final (default: 1.0)
  -s EPS        Softening (default: 0.01)
  -o FILE       Archivo de salida CSV
  --init TYPE   Tipo de init: plummer|uniform|twobody
  --validate    Correr suite de validacion
  --seed SEED   Semilla aleatoria (default: 42)
```

## Compilacion

```bash
make          # compila nbody con gcc -O2
make clean    # limpia binarios
```

## Siguiente paso: Semana 2

La version secuencial esta lista como base para la semana 2, donde se implementara Barnes-Hut (octree + calculo de fuerza O(N log N)).
