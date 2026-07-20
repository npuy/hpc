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
 * Recorre el árbol en post-orden llenando mass y cm de cada nodo: las hojas
 * agregan sus partículas y los nodos internos combinan la masa y el centro de
 * masa ponderado de sus hijos. Debe llamarse tras octree_build y antes de
 * calcular fuerzas.
 */
void octree_compute_mass(Octree *t, const Particle *p);

/* Libera el pool de nodos, el arreglo next[] y el propio Octree. */
void octree_free(Octree *t);

#endif
