# Semana 1 — Plan Detallado: Version Secuencial

## Estructura de archivos objetivo

```
src/
  main.c              -- programa principal, loop de simulacion
  particle.h / .c     -- struct Particle, inicializacion, I/O
  force.h / .c        -- calculo de fuerza directa O(N^2)
  leapfrog.h / .c     -- integrador leapfrog (kick-drift-kick)
  morton.h / .c       -- codificacion Morton 3D
  vec3.h              -- macros/inline para operaciones vectoriales
  validation.h / .c   -- tests de conservacion y correctitud
Makefile
```

## Dia 1-2: Particulas + Infraestructura

### Estudiante A

#### 1. particle.h/.c — Estructura base

```c
typedef struct {
    double pos[3];    // x, y, z
    double vel[3];    // vx, vy, vz
    double acc[3];    // ax, ay, az
    double mass;
    uint64_t morton;  // se llena despues
} Particle;
```

#### 2. Generadores de condiciones iniciales

- `init_two_body(Particle *p, int *n)` — Dos cuerpos en orbita circular (caso mas simple para validar).
- `init_plummer(Particle *p, int n, unsigned seed)` — Distribucion de Plummer (modelo astrofisico estandar para N cuerpos). Generar masa acumulada uniforme, derivar radio, asignar velocidad con rejection sampling.
- `init_uniform_cube(Particle *p, int n, unsigned seed)` — Distribucion uniforme en cubo (util para testing).

#### 3. I/O

- `write_particles_csv(const char *fname, Particle *p, int n, double t)` — Salida CSV para analisis/plotting.
- `read_particles_csv(const char *fname, Particle *p, int *n)` — Lectura para restart.

#### 4. vec3.h — Macros inline

```c
#define VEC3_SUB(r, a, b) do { (r)[0]=(a)[0]-(b)[0]; (r)[1]=(a)[1]-(b)[1]; (r)[2]=(a)[2]-(b)[2]; } while(0)
#define VEC3_ADD(r, a, b) do { (r)[0]=(a)[0]+(b)[0]; (r)[1]=(a)[1]+(b)[1]; (r)[2]=(a)[2]+(b)[2]; } while(0)
#define VEC3_SCALE(r, a, s) do { (r)[0]=(a)[0]*(s); (r)[1]=(a)[1]*(s); (r)[2]=(a)[2]*(s); } while(0)
#define VEC3_DOT(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define VEC3_NORM(a) sqrt(VEC3_DOT(a, a))
```

#### 5. Makefile

- Target `nbody` con `gcc -O2 -Wall -lm`
- Target `nbody_mpi` (para futuro, con `mpicc`)
- Target `clean`

### Estudiante B

#### 1. force.h/.c — Fuerza directa O(N^2)

```c
void compute_forces_direct(Particle *p, int n, double softening);
```

Para cada par (i,j):

```
r = pos[j] - pos[i]
dist = sqrt(r.r + eps^2)       // eps = softening
f = G * m[i] * m[j] / dist^3
acc[i] += f * r / m[i]
```

- **Softening** (`eps ~ 0.01`): evita singularidades cuando dos particulas estan muy cerca.
- Usar simetria de Newton (F_ij = -F_ji) para reducir calculos a N(N-1)/2.
- **G = 1** (unidades naturales).

#### 2. Funciones de energia (para validacion)

```c
double kinetic_energy(Particle *p, int n);
double potential_energy(Particle *p, int n, double softening);
```

## Dia 3-4: Leapfrog + Morton

### Estudiante A

#### 1. leapfrog.h/.c — Integrador Kick-Drift-Kick (KDK)

```c
void leapfrog_kick(Particle *p, int n, double dt_half);
void leapfrog_drift(Particle *p, int n, double dt);
```

- **Kick**: `vel += acc * dt/2`
- **Drift**: `pos += vel * dt`
- El loop principal hace: kick(dt/2) -> drift(dt) -> compute_forces -> kick(dt/2)
- En el paso 0: calcular fuerzas iniciales antes del primer kick.

KDK permite tener posiciones y velocidades sincronizadas al final de cada paso, lo que facilita calcular energia correctamente.

#### 2. morton.h/.c — Codificacion Morton 3D

```c
uint64_t morton_encode(double x, double y, double z,
                       double min[3], double max[3], int bits);
void morton_sort(Particle *p, int n, double min[3], double max[3]);
```

- **Normalizar** posiciones al rango [0, 2^bits) usando el bounding box.
- **Interleaving**: expandir cada coordenada a 21 bits (63 bits total en uint64), intercalar bits x,y,z.
- **Ordenar** particulas por clave Morton (qsort).
- Usar `bits = 21` para maxima resolucion en 64 bits (21 * 3 = 63).

Funcion auxiliar de bit-interleaving:

```c
static uint64_t expand_bits(uint64_t v) {
    v &= 0x1fffff;  // 21 bits
    v = (v | v << 32) & 0x1f00000000ffff;
    v = (v | v << 16) & 0x1f0000ff0000ff;
    v = (v | v << 8)  & 0x100f00f00f00f00f;
    v = (v | v << 4)  & 0x10c30c30c30c30c3;
    v = (v | v << 2)  & 0x1249249249249249;
    return v;
}
```

### Estudiante B

#### validation.h/.c — Suite de validacion

**Test 1: Orbita circular de dos cuerpos**

- Inicializar dos masas iguales en orbita circular.
- Correr 1 periodo orbital completo.
- Verificar que la posicion final coincide con la inicial (error < 1% con dt suficientemente pequeno).

**Test 2: Conservacion de energia**

- Con N=100 particulas (Plummer), correr 100 pasos.
- Verificar que |E(t) - E(0)| / |E(0)| < 1e-3 (leapfrog es simplectico, conserva energia bien).

**Test 3: Conservacion de momento lineal**

- Verificar que el momento total del sistema se conserva (deberia ser exacto a nivel de redondeo).

**Test 4: Morton ordering**

- Crear 8 puntos en las esquinas de un cubo unitario.
- Verificar que las claves Morton estan en el orden esperado (Z-order).
- Verificar que puntos espacialmente cercanos tienen claves cercanas.

## Dia 5: Integracion + Main Loop

### Ambos

#### main.c — Programa principal

```c
int main(int argc, char *argv[]) {
    int N = 1000;           // default
    double dt = 0.001;
    double t_end = 1.0;
    double softening = 0.01;

    // Parsear argumentos: -n N -dt DT -t T_END -o output

    Particle *particles = malloc(N * sizeof(Particle));
    init_plummer(particles, N, 42);

    // Fuerzas iniciales
    compute_forces_direct(particles, N, softening);

    double E0 = kinetic_energy(particles, N) + potential_energy(particles, N, softening);

    for (double t = 0; t < t_end; t += dt) {
        leapfrog_kick(particles, N, dt / 2);
        leapfrog_drift(particles, N, dt);
        compute_forces_direct(particles, N, softening);
        leapfrog_kick(particles, N, dt / 2);

        // Cada M pasos: escribir snapshot, verificar energia
    }

    double Ef = kinetic_energy(particles, N) + potential_energy(particles, N, softening);
    printf("Energy drift: %.6e\n", fabs(Ef - E0) / fabs(E0));

    free(particles);
}
```

#### Argumentos CLI

| Flag       | Descripcion                                          |
| ---------- | ---------------------------------------------------- |
| `-n`       | numero de particulas                                 |
| `-dt`      | paso temporal                                        |
| `-t`       | tiempo final                                         |
| `-s`       | softening                                            |
| `-o`       | archivo de salida                                    |
| `--init`   | tipo de inicializacion (plummer/uniform/twobody)     |
| `--validate` | correr suite de validacion en vez de simulacion    |

## Criterios de aceptacion

| Test                                       | Criterio                       |
| ------------------------------------------ | ------------------------------ |
| Dos cuerpos completan orbita               | Error posicion < 1%            |
| Conservacion energia (N=1000, 1000 pasos)  | Drift < 0.1%                   |
| Conservacion momento lineal                | Drift < 1e-12                  |
| Morton: 8 esquinas del cubo                | Orden correcto                 |
| Morton sort + verificacion                 | Particulas ordenadas           |
| Simulacion N=10000 completa sin crash      | ---                            |
| Tiempo razonable N=10000 O(N^2)            | < 30s por paso en laptop       |

## Parametros recomendados

| Parametro  | Valor dev | Valor test |
| ---------- | --------- | ---------- |
| N          | 100-1000  | 10000      |
| dt         | 0.001     | 0.0005     |
| softening  | 0.01      | 0.01       |
| t_end      | 1.0       | 10.0       |
| G          | 1.0       | 1.0        |

## Notas tecnicas

- **Unidades naturales**: G=1, masas ~1/N (masa total = 1), posiciones ~1.
- **Softening**: sin el, particulas cercanas generan fuerzas infinitas y la simulacion explota. eps = 0.01 es un buen default.
- **Leapfrog KDK vs DKD**: KDK (kick-drift-kick) permite tener posiciones y velocidades sincronizadas al final de cada paso, lo que facilita calcular energia correctamente.
- **Compilacion local (Mac)**: `gcc -O2 -Wall -lm -o nbody src/*.c` (sin MPI por ahora).
- **Compilacion en fing**: `mpicc -O2 -Wall -lm -o nbody src/*.c` (cuando se agregue MPI).
