# Spec: ORB (bisección recursiva ortogonal)

Semana 5. Reemplaza la descomposición por rangos de clave Morton.

## Por qué

La semana 4 midió que el LET no comprime sobre Plummer: `n_ghost` crecía como
N^0,92, la misma asintótica que la réplica. La causa medida es que **un tramo
contiguo de claves Morton no es una región compacta** — el AABB del dominio de un
proceso llegó a cubrir el 44,6% del espacio, con lo que ningún nodo satisfacía el
criterio de exportación.

El problema de fondo: en Plummer el pico de densidad cae justo donde se tocan los
octantes de nivel 1, así que las claves del core se reparten por todo el espacio
de claves y cualquier corte por percentil le da a cada proceso una mezcla de core
y halo. Los dominios se interpenetran.

ORB corta el espacio por planos: cada proceso queda con una **caja**, los dominios
son disjuntos y contiguos, y la caja de un proceso periférico no contiene el core.
Recién ahí el criterio del LET puede aceptar resúmenes.

## Estructura

```c
typedef struct {
    double min[3], max[3];  /* caja de este nodo */
    int    first, nprocs;   /* tramo de ranks [first, first+nprocs) */
    int    axis;            /* eje de corte; -1 si es hoja */
    double split;           /* coordenada del corte */
    int    left, right;     /* indices de hijos; -1 si es hoja */
} OrbNode;
```

El árbol tiene `2P-1` nodos y está **replicado en todos los procesos**: son unos
cientos de bytes, así que cada proceso conoce la caja de todos los demás sin
comunicación. Es justo lo que el LET necesita.

## Algoritmo

```
orb_build(caja, idx[], [first, first+nprocs)):
    si nprocs == 1: hoja, first es el dueño
    axis = dimensión más larga de la caja
    w_total = Allreduce(Σ peso de idx[])
    nL = nprocs / 2
    objetivo = w_total * nL / nprocs
    bisecar c en [caja.min[axis], caja.max[axis]], 50 veces:
        below = Allreduce(Σ peso de idx[] con pos[axis] < c)
        si below < objetivo: lo = c   si no: hi = c
    particionar idx[] in-place por pos[axis] < split
    recurrir en las dos mitades
```

### Decisiones

- **Eje = dimensión más larga.** Evita cajas laminares, que tienen mucha
  superficie y por lo tanto generan mucho LET.
- **Se bisecta el trabajo, no el conteo.** El campo `work` (semana 4) ya estaba
  validado, así que ORB hereda el balance por costo sin código nuevo. Con `work`
  global nulo (primer paso) cae al conteo.
- **50 iteraciones fijas, no tolerancia relativa.** Todos los procesos deben hacer
  exactamente las mismas llamadas colectivas; un criterio de corte que dependa de
  datos locales rompería el lockstep.
- **P que no es potencia de 2** sale con reparto proporcional, sin parches.
  Verificado con P=3.

### La pertenencia se lleva por el camino, no por geometría

El arreglo de índices locales se particiona in-place en cada nivel con el **mismo
predicado `<`** que usa `domain_owner_pos`. No hay test de contención geométrica.

Esto vuelve imposible por construcción el riesgo señalado en el plan: si el conteo
de la bisección y la asignación de dueño usaran comparaciones distintas, una
partícula justo sobre el plano se contaría de un lado y se enviaría al otro, y se
perdería. Llevarlo por el camino recorrido lo hace estructuralmente imposible en
vez de dependiente de que dos comparaciones coincidan.

Ventaja secundaria: cada nivel cuesta O(n_idx) en total, no O(n_local) por nodo.

## API

```c
void domain_partition_orb(Domain *d, const Particle *p, int n_local,
                          int use_work, MPI_Comm comm);   /* colectivo */
int  domain_owner_pos(const Domain *d, const double pos[3]);  /* O(log P) */
void domain_box(const Domain *d, int rank, double min[3], double max[3]);
```

`domain_owner_pos` **no hace test de contención**: cualquier posición, aun fuera
de la caja global, obtiene un dueño válido. Es lo que hace seguro el caso de una
partícula que se escapó del box congelado entre reparticiones.

## Morton cambia de rol

Ya no determina la propiedad. Queda **solo para localidad de caché** al ordenar el
arreglo local antes de construir el árbol.

| Invariante | Estado |
|---|---|
| `replicate_particles` devuelve un arreglo globalmente ordenado por Morton | **Deja de valer con ORB.** Los tests que lo usan emparejan por `id`. |
| Arreglo local ordenado por clave para la búsqueda binaria de `domain_partition` | Solo aplica al camino Morton |
| `Domain::gmin/gmax` congelado | **Sigue haciendo falta**, ahora porque `octree_build_box` exige la misma caja raíz en todos los procesos |
| `morton_encode` acota las coordenadas (bug de la semana 3) | Sigue vigente: las claves se siguen calculando |

## Resultados

| Métrica | Morton | ORB |
|---|---|---|
| Exponente de `n_ghost` ~ N^k | 0,97 | **0,66** |
| `tree_time` de P=1 a P=8 | 4,8× | **1,6×** |
| Fantasmas (N=100K, P=4) | 61.056 | **16.933** |
| Brecha Plummer/uniform | 2,77× | **1,83×** |

El exponente 0,66 está muy cerca del 2/3 de un efecto de superficie, que es lo que
un LET debe costar.

## Comparación conservada

`--decomp orb|morton` mantiene los dos caminos vivos. Es lo que permitió construir
todas las tablas comparativas del reporte, y el test 12 valida ambos.
