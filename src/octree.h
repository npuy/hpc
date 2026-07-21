#ifndef OCTREE_H
#define OCTREE_H

/*
 * octree.h — Octree espacial para el algoritmo de Barnes-Hut.
 *
 * Subdivide recursivamente el dominio cúbico en 8 octantes hasta que cada
 * hoja contiene a lo sumo una partícula. Cada nodo almacena la masa total y
 * el centro de masa de su subárbol, lo que permite aproximar la fuerza de
 * grupos lejanos por un único cuerpo equivalente (ver barnes_hut.h).
 *
 * Se usa un pool de nodos (arreglo dinámico) con índices int en lugar de
 * punteros: reduce llamadas a malloc, mejora la localidad de caché y es
 * directamente serializable para el intercambio LET por MPI (semanas 3-4).
 * IMPORTANTE: el pool puede reubicarse por realloc; nunca cachear punteros
 * a nodes[] a través de una inserción — usar siempre índices.
 */

#include "particle.h"

/*
 * Nodo del octree. Un nodo puede estar en tres estados:
 *   - hoja vacía:     is_leaf=1, particle=-1
 *   - hoja ocupada:   is_leaf=1, particle>=0 (cabeza de una cadena en next[])
 *   - nodo interno:   is_leaf=0, children[] apuntando a los 8 octantes
 * La cadena via next[] solo aparece para partículas coincidentes o al llegar
 * a la profundidad máxima; en el caso normal cada hoja tiene una partícula.
 */
typedef struct {
    double min[3], max[3];   /* bounding box cúbico de la celda */
    double mass;             /* masa total del subárbol */
    double cm[3];            /* centro de masa del subárbol */
    int    children[8];      /* índices a los 8 hijos; -1 si no subdividido */
    int    particle;         /* índice de partícula (cabeza de cadena) o -1 */
    int    is_leaf;          /* 1 si es hoja, 0 si es nodo interno */
} OctreeNode;

typedef struct {
    OctreeNode *nodes;       /* pool de nodos */
    int         count;       /* nodos usados */
    int         capacity;    /* nodos asignados */
    int        *next;        /* cadena de partículas por hoja (tam n); -1 = fin */
    int         root;        /* índice de la raíz */
    double      size;        /* lado del cubo raíz (para el criterio s/d) */
} Octree;

/*
 * Calcula un bounding box cúbico que contiene a las n partículas, centrado en
 * el centro del AABB y con un margen del 1% para que ninguna partícula caiga
 * exactamente sobre max[k]. Escribe las esquinas en min[3] y max[3].
 */
void octree_bounds(const Particle *p, int n, double min[3], double max[3]);

/*
 * Construye el octree insertando las n partículas una a una en el cubo global.
 * Reserva y retorna un Octree nuevo (el caller debe liberar con octree_free).
 * No calcula masas/centros de masa; llamar a octree_compute_mass después.
 */
Octree *octree_build(const Particle *p, int n);

/*
 * Igual que octree_build pero con el cubo raíz dado EXPLÍCITAMENTE en vez de
 * derivarlo de las partículas recibidas.
 *
 * Es la variante que exige el LET. octree_build calcula el cubo a partir de las
 * partículas que le pasan: con replicación todos los procesos pasaban el mismo
 * arreglo global y obtenían el mismo cubo, pero con LET cada proceso construye
 * sobre sus N/P locales y obtendría un cubo distinto. Las celdas no alinearían
 * entre procesos y las aproximaciones importadas no serían comparables con las
 * locales. Todos los procesos deben pasar el mismo cubo (Domain::gmin/gmax).
 *
 * Precondición: las partículas deberían caer dentro del cubo. Una que se salga
 * no rompe la estructura (siempre hay un octante y la recursión termina en
 * MAX_DEPTH), pero queda en una celda que no la contiene, y entonces el lado s
 * del nodo subestima la extensión real de su contenido y el criterio de
 * apertura puede quedarse corto. Es el mismo supuesto que ya hace morton_encode
 * al acotar, y lo que lo mantiene válido es repartir cada K pasos.
 */
Octree *octree_build_box(const Particle *p, int n,
                         const double min[3], const double max[3]);

/*
 * Inserta en un árbol ya construido las partículas p[n_old..n_new), donde p es
 * el arreglo extendido (las primeras n_old deben ser las mismas y estar en el
 * mismo orden con que se construyó t).
 *
 * Existe para el LET: permite agregar los fantasmas importados al árbol local
 * en vez de reconstruir un árbol fusionado desde cero. La caja raíz no cambia,
 * así que la estructura ya calculada sigue siendo válida y solo se subdividen
 * las hojas donde caen los nuevos puntos. Ahorra una construcción completa por
 * paso, que es justamente la fase que la semana 4 busca abaratar.
 *
 * Hay que volver a llamar a octree_compute_mass después: las masas y centros de
 * masa del árbol quedan desactualizados al agregar puntos.
 */
void octree_insert_particles(Octree *t, const Particle *p, int n_old, int n_new);

/*
 * Recorre el árbol en post-orden llenando mass y cm de cada nodo: las hojas
 * agregan sus partículas y los nodos internos combinan la masa y el centro de
 * masa ponderado de sus hijos. Debe llamarse tras octree_build y antes de
 * calcular fuerzas.
 */
void octree_compute_mass(Octree *t, const Particle *p);

/* Libera el pool de nodos, el arreglo next[] y el propio Octree. */
void octree_free(Octree *t);

#endif
