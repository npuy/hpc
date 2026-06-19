# Módulo: particle

## Descripción general

Define la estructura de datos fundamental `Particle` y provee funciones para inicializar distribuciones de partículas y para persistirlas en formato CSV. Es el módulo central del que dependen todos los demás.

## Estructura de datos

```c
typedef struct {
    double pos[3];    // Posición en espacio 3D (x, y, z)
    double vel[3];    // Velocidad (vx, vy, vz)
    double acc[3];    // Aceleración (ax, ay, az)
    double mass;      // Masa de la partícula
    uint64_t morton;  // Clave Morton para ordenamiento Z-order
} Particle;
```

**Layout de memoria**: 9 doubles (72 bytes) + 1 uint64_t (8 bytes) = 80 bytes por partícula. El campo `morton` se calcula externamente por el módulo `morton`.

## Interfaz pública

### Funciones de inicialización

| Función | Descripción |
|---|---|
| `init_two_body(p, &n)` | Configura un sistema de 2 cuerpos de masa igual (`M/2` cada uno) en órbita circular simétrica respecto al origen. Separación unitaria, velocidades calculadas analíticamente para órbita circular con G=1. Escribe `n=2`. |
| `init_plummer(p, n, seed)` | Genera `n` partículas según el modelo de Plummer (esfera isotérmica truncada). Las posiciones se muestrean con la CDF inversa y las velocidades con rejection sampling (Aarseth, Hénon, Wielen 1974). Aplica corrección de centro de masa para que el sistema quede centrado en el origen con momento total nulo. Masa total = 1, radio de escala `a = 1`. |
| `init_uniform_cube(p, n, seed)` | Genera `n` partículas con posiciones uniformes en el cubo `[-0.5, 0.5]³`. Velocidades iniciales nulas. Masa total = 1 repartida equitativamente. |

### Funciones de E/S

| Función | Descripción |
|---|---|
| `write_particles_csv(fname, p, n, t)` | Escribe el estado de `n` partículas al archivo `fname` en formato CSV con header `t,mass,x,y,z,vx,vy,vz`. Precisión: 12 dígitos significativos en notación científica. |
| `read_particles_csv(fname, p, &n)` | Lee partículas desde un CSV con el formato anterior. Salta la línea de header. Retorna `0` en éxito, `-1` si no puede abrir el archivo. Escribe la cantidad leída en `*n`. |

## Formato CSV

```
t,mass,x,y,z,vx,vy,vz
1.000000e+00,5.000000e-01,1.234567e-01,...
```

Cada fila contiene el tiempo del snapshot, la masa, 3 coordenadas de posición y 3 de velocidad.

## Dependencias

- `vec3.h` — usado internamente pero no expuesto en la interfaz.
- `<stdint.h>` — para `uint64_t` en el campo `morton`.

## Limitaciones

- `read_particles_csv` asume que el caller ya alocó un buffer `Particle*` suficientemente grande; no hay protección contra overflow.
- Las funciones de inicialización usan `rand()` (no thread-safe); en la versión MPI/OpenMP se deberá reemplazar por un PRNG independiente por hilo.
