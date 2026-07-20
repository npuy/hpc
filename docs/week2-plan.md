# Semana 2 — Plan Detallado: Barnes-Hut Secuencial

## Objetivo

Reemplazar el cálculo de fuerzas directo O(N²) por el algoritmo de **Barnes-Hut**
O(N log N): construir un **octree** que agrupa partículas espacialmente y aproximar
la fuerza de grupos lejanos por su centro de masa. Todo sigue siendo secuencial;
la meta es correctitud y escalabilidad algorítmica, no paralelismo (eso es semana 3).

Corresponde a las **etapas 2 y 4** del plan general:

| Etapa | Objetivo         | Validación                    |
| ----- | ---------------- | ----------------------------- |
| 4     | Octree           | Masa total, centros de masa   |
| 2     | Barnes-Hut local | Error respecto a O(N²)        |

**Precondición (semana 1, ya cumplida):** versión secuencial validada con fuerza
directa, leapfrog KDK, Morton sort y suite de validación (4/4 tests PASS).

## Estado de partida (qué existe)

```
src/
  vec3.h              -- macros vectoriales 3D           [reutilizar]
  particle.h / .c     -- struct Particle (con .morton), init, I/O   [reutilizar]
  force.h / .c        -- compute_forces_direct O(N^2), energias      [oraculo de referencia]
  leapfrog.h / .c     -- integrador KDK                  [reutilizar sin cambios]
  morton.h / .c       -- morton_encode, morton_sort      [base del octree]
  validation.h / .c   -- 4 tests secuenciales            [extender]
  main.c              -- loop de simulacion              [agregar seleccion de metodo]
```

`compute_forces_direct` se conserva intacto: es el **oráculo** contra el cual se mide
el error de Barnes-Hut.

## Estructura de archivos objetivo (nuevos)

```
src/
  octree.h / .c       -- construccion del octree, centros de masa, bounding box
  barnes_hut.h / .c   -- calculo de fuerza por traversal con criterio de apertura (MAC)
docs/specs/
  octree.md           -- spec del modulo octree
  barnes_hut.md       -- spec del modulo barnes_hut
```

## Estructuras de datos

### OctreeNode (array-pool, no punteros)

Se usa un **pool de nodos** (arreglo dinámico `OctreeNode *nodes`) e índices `int`
en lugar de punteros. Ventajas: menos `malloc`/`free`, mejor localidad de caché, y
—clave para semanas 3-4— serializable directamente para el intercambio LET por MPI.

```c
typedef struct {
    double min[3], max[3];   /* bounding box de la celda (cubo) */
    double mass;             /* masa total contenida en el subarbol */
    double cm[3];            /* centro de masa del subarbol */
    int    children[8];      /* indices a nodos hijos; -1 si vacio */
    int    particle;         /* indice de particula si is_leaf; -1 si interno/vacio */
    int    is_leaf;          /* 1 si es hoja (0 o 1 particula) */
} OctreeNode;

typedef struct {
    OctreeNode *nodes;       /* pool de nodos */
    int         count;       /* nodos usados */
    int         capacity;    /* nodos asignados */
    int         root;        /* indice de la raiz (normalmente 0) */
    double      size;        /* lado del cubo raiz (para el criterio s/d) */
} Octree;
```

**Convención de octantes** (índice del hijo por bits del signo respecto al centro):

```
child = (x >= cx ? 1 : 0) | (y >= cy ? 2 : 0) | (z >= cz ? 4 : 0)
```

## Reparto de trabajo

- **Estudiante A → `octree.h/.c`** (construcción, bounding box, centros de masa).
- **Estudiante B → `barnes_hut.h/.c`** (traversal, criterio de apertura, integración
  con `main.c`) y **extensión de la suite de validación**.

La interfaz entre ambos es el struct `Octree` y sus funciones públicas; acordarla el
Día 1 para poder trabajar en paralelo.

---

## Día 1-2: Construcción del octree

### Estudiante A

#### 1. Bounding box global (cubo)

```c
void octree_bounds(const Particle *p, int n, double min[3], double max[3]);
```

- Calcular el AABB de todas las partículas.
- **Convertirlo en cubo** (lado = máx de las 3 extensiones) y agregar un margen
  del ~1% para que ninguna partícula caiga exactamente sobre `max[k]`. Un cubo
  simétrico garantiza que la subdivisión en octantes sea uniforme.

#### 2. Construcción top-down recursiva

```c
Octree *octree_build(const Particle *p, int n);
void    octree_free(Octree *t);
```

Algoritmo recomendado (inserción incremental con subdivisión perezosa):

```
crear raiz con el bounding box cubico global
para cada particula i:
    octree_insert(t, root, p, i)

octree_insert(nodo, particula i):
    si nodo vacio            -> guardar i como hoja
    si nodo es hoja (tiene j)-> subdividir: crear 8 hijos,
                                reinsertar j y luego i en el octante que corresponda
    si nodo interno          -> descender al octante de i e insertar recursivamente
```

- `octree_alloc_node(t)` crece el pool con realloc (duplicando capacidad).
  **Cuidado:** guardar índices, no punteros, porque `realloc` puede mover `nodes`.
- **Caso degenerado:** dos partículas en la misma posición (o a distancia menor que
  la celda mínima) producen recursión infinita. Cortar la subdivisión a una
  profundidad máxima (`MAX_DEPTH ~ 21`, límite de Morton) y permitir hojas con
  bucket de >1 partícula, o desplazar con el softening. Documentar la decisión.

> **Alternativa (opcional, más rápida):** construir bottom-up aprovechando que las
> partículas ya vienen ordenadas por clave Morton (semana 1). Los prefijos comunes
> de las claves definen la jerarquía de celdas. Es más eficiente pero más complejo;
> dejar como mejora si sobra tiempo. La versión por inserción es suficiente y más
> fácil de depurar.

#### 3. Centros de masa (post-orden)

```c
void octree_compute_mass(Octree *t);   /* llena mass[] y cm[] de todos los nodos */
```

Recorrido post-orden desde la raíz:

```
hoja con particula j -> mass = m_j; cm = pos_j
nodo interno         -> mass = Σ mass_hijo
                        cm   = (Σ mass_hijo * cm_hijo) / mass
```

Se calcula **una vez por paso**, después de construir el árbol y antes de la fuerza.

### Entregable Día 2 (A)

- `octree_build` produce un árbol sin crashes para N=1000 Plummer.
- `octree_compute_mass` con **masa de la raíz == masa total del sistema** (a
  redondeo) y **cm de la raíz == centro de masa del sistema**.

---

## Día 3-4: Cálculo de fuerza Barnes-Hut

### Estudiante B

#### 1. Criterio de apertura (MAC) y traversal

```c
void compute_forces_bh(const Octree *t, Particle *p, int n,
                       double theta, double softening);
```

Para cada partícula `i`, recorrer el árbol desde la raíz:

```
fuerza_sobre(i, nodo):
    si nodo vacio: retornar
    si nodo es hoja:
        si la particula de la hoja != i: acumular fuerza directa (2 cuerpos)
    si nodo interno:
        s = lado de la celda del nodo
        d = distancia de i al cm del nodo
        si (s / d) < theta:                     # nodo suficientemente lejano
            acumular fuerza usando (mass, cm) del nodo como un solo cuerpo
        si no:
            recursar en los 8 hijos
```

- Reutilizar la **misma fórmula de fuerza con softening** que `force.c` para que el
  límite `theta → 0` coincida numéricamente con el método directo:
  `f = G · m_i · m_nodo / (d² + ε²)^{3/2} · r̂`.
- Resetear `acc[]` a cero al inicio (igual que `compute_forces_direct`).
- **θ (theta):** parámetro de apertura. `θ = 0` ⇒ exacto (equivale a O(N²));
  `θ ≈ 0.5` es el estándar (buen balance error/velocidad); `θ` mayor ⇒ más rápido y
  menos preciso. Precomputar `s` por nodo o derivarlo de `t->size / 2^depth`.

#### 2. Energía potencial con Barnes-Hut (opcional pero recomendado)

Para monitorear conservación de energía sin volver a O(N²) en N grande, una variante
`potential_energy_bh` que use el mismo traversal. Para validación con N chico, basta
la `potential_energy` directa existente.

### Estudiante B — Extensión de la suite de validación

Agregar a `validation.c` (invocables con `--validate`):

**Test 5: Conservación de masa del árbol**
- Construir octree para N=1000 Plummer.
- `|masa_raiz − masa_total| / masa_total < 1e-12`.

**Test 6: Centro de masa del árbol**
- `‖cm_raiz − cm_directo‖ < 1e-10` respecto al CM calculado por suma directa.

**Test 7: Error de fuerza Barnes-Hut vs directo**
- N=1000, calcular `acc` con ambos métodos.
- Error relativo L2 por partícula:
  `err = ‖a_bh − a_dir‖ / ‖a_dir‖`, promediado y máximo.
- Criterio: `err_medio < 1e-2` para `θ = 0.5`.

**Test 8: Convergencia en θ**
- Barrer `θ ∈ {0.8, 0.5, 0.2, 0.0}` y verificar que el error **decrece
  monótonamente** y que `θ = 0` reproduce el método directo (`err < 1e-9`).

---

## Día 5: Integración en el loop y benchmarking

### Ambos

#### 1. `main.c` — selección de método

Nuevos argumentos CLI:

| Flag         | Descripción                                    | Default  |
| ------------ | ---------------------------------------------- | -------- |
| `--method`   | `direct` \| `bh`                               | `bh`     |
| `--theta`    | ángulo de apertura Barnes-Hut                  | `0.5`    |

Loop de simulación con Barnes-Hut (el árbol **se reconstruye cada paso** porque las
partículas se mueven):

```c
for (double t = 0; t < t_end; t += dt) {
    leapfrog_kick(particles, N, dt / 2);
    leapfrog_drift(particles, N, dt);

    if (method == BH) {
        Octree *tree = octree_build(particles, N);
        octree_compute_mass(tree);
        compute_forces_bh(tree, particles, N, theta, softening);
        octree_free(tree);
    } else {
        compute_forces_direct(particles, N, softening);
    }

    leapfrog_kick(particles, N, dt / 2);
    // cada M pasos: snapshot + verificar energia
}
```

- Opcional: reordenar por Morton cada K pasos para mantener localidad de caché a
  medida que las partículas migran (`morton_sort` ya existe).
- Instrumentar por separado el **tiempo de construcción del árbol** vs **tiempo de
  fuerza** (útil para el informe y para semana 3).

#### 2. Benchmark de escalabilidad O(N log N)

Correr `direct` vs `bh` (θ=0.5) para N ∈ {1K, 5K, 10K, 50K, 100K}, 10 pasos:

- Verificar que el tiempo de `bh` crece ~N log N mientras `direct` crece ~N².
- Reportar el **N de cruce** donde Barnes-Hut se vuelve más rápido.
- Registrar drift de energía de `bh` (θ=0.5) sobre 1000 pasos (N=1000): debe seguir
  siendo pequeño (la aproximación introduce ruido, pero leapfrog lo mantiene acotado).

---

## Criterios de aceptación

| Test / Criterio                                          | Umbral                       |
| -------------------------------------------------------- | ---------------------------- |
| Masa de la raíz == masa total                            | error < 1e-12                |
| CM de la raíz == CM directo                              | ‖·‖ < 1e-10                  |
| Error fuerza BH vs directo (N=1000, θ=0.5)               | err medio < 1%               |
| Convergencia θ→0 reproduce método directo                | err < 1e-9 en θ=0            |
| Error decrece monótonamente al bajar θ                   | monótono                     |
| Conservación de energía con BH (N=1000, 1000 pasos, θ=0.5) | drift < 1%                 |
| Escalabilidad: BH más rápido que directo para N≥50K      | tiempo_bh < tiempo_direct    |
| Simulación N=100K con BH completa sin crash              | ---                          |
| Sin fugas de memoria (árbol reconstruido cada paso)      | `valgrind` limpio            |

## Parámetros recomendados

| Parámetro   | Valor dev | Valor test |
| ----------- | --------- | ---------- |
| N           | 1K        | 100K       |
| theta       | 0.5       | 0.5        |
| dt          | 0.001     | 0.0005     |
| softening   | 0.01      | 0.01       |
| t_end       | 1.0       | 1.0        |

## Notas técnicas

- **Por qué octree y no kd-tree:** el octree (8 hijos, subdivisión regular en 3D)
  encaja naturalmente con la codificación Morton ya implementada y con la
  descomposición de dominio por rangos de clave que se usará en MPI (semanas 3-4).
- **Reconstrucción por paso:** en esta etapa el árbol se rearma completo cada paso.
  Es lo correcto y simple; optimizaciones (árbol incremental, reuso) quedan fuera de
  alcance. Medir el costo relativo de construcción vs fuerza.
- **Array-pool vs punteros:** usar índices `int` en `children[8]` (no punteros)
  evita fugas, mejora la localidad y prepara la serialización para LET/MPI. Recordar
  no cachear punteros a `nodes[]` a través de un `realloc`.
- **Estabilidad numérica:** usar exactamente la misma fórmula con softening que
  `force.c`; así θ=0 debe coincidir con el método directo salvo orden de sumas.
- **Casos degenerados:** partículas coincidentes o cúmulos muy densos ⇒ limitar
  `MAX_DEPTH` y/o permitir buckets en hojas para evitar recursión infinita.
- **Compilación local (Mac):** `gcc -O2 -Wall -lm -o nbody src/*.c` (sigue sin MPI).
- **Actualizar Makefile:** agregar `octree.c` y `barnes_hut.c` a la lista de fuentes.

## Riesgos de la semana

| Riesgo                          | Mitigación                                       |
| ------------------------------- | ------------------------------------------------ |
| Error Barnes-Hut (etapa)        | Comparar acc contra `compute_forces_direct`      |
| Recursión infinita (degenerado) | `MAX_DEPTH` + buckets en hoja                     |
| Fuga de memoria por paso        | `octree_free` + chequeo con `valgrind`            |
| Puntero inválido tras `realloc` | Usar índices, nunca cachear `&nodes[i]`           |
| θ mal elegido (lento o impreciso)| Barrido de θ en test de convergencia             |

## Entregables de la semana

1. `src/octree.h/.c` y `src/barnes_hut.h/.c` compilando y sin fugas.
2. `main.c` con `--method` y `--theta`; Barnes-Hut integrado en el loop.
3. Tests 5-8 agregados a la suite (`--validate` pasa 8/8).
4. `docs/specs/octree.md` y `docs/specs/barnes_hut.md`.
5. `docs/week2-report.md` con tabla de validación y benchmark direct vs BH.
