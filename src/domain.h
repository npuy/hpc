#ifndef DOMAIN_H
#define DOMAIN_H

/*
 * domain.h — Descomposición del dominio por rangos de clave Morton.
 *
 * El espacio de claves [0, 2^63) se corta en nprocs tramos contiguos, uno por
 * proceso. Como la curva Z preserva localidad, un tramo contiguo de claves es
 * una región compacta del espacio 3D: cada proceso queda con un vecindario
 * espacial, que es lo que hace baratos el LET y la migración.
 *
 * Los cortes NO se hacen en partes iguales del espacio de claves. Una
 * distribución de Plummer concentra casi todas las partículas en unas pocas
 * celdas centrales, así que repartir por volumen dejaría a un proceso con casi
 * todo el trabajo. Se reparte por CONTEO de partículas, buscando los P-1
 * splitters que dejan N/P partículas en cada tramo.
 *
 * El algoritmo es una bisección sobre el espacio de claves: en cada iteración
 * se cuenta cuántas partículas globales caen por debajo de un candidato
 * (MPI_Allreduce) y se ajusta el intervalo de búsqueda. Los P-1 splitters se
 * bisectan SIMULTÁNEAMENTE, de modo que el costo total es de ~63 colectivos y
 * no de 63·(P-1). No requiere ordenar globalmente ni juntar las claves en un
 * proceso, y es exactamente el procedimiento que reusará el rebalanceo dinámico
 * de la semana 4: allí basta con volver a llamarlo cada K pasos.
 */

#ifdef USE_MPI

#include <mpi.h>
#include <stdint.h>
#include "particle.h"

/*
 * Nodo del árbol ORB (bisección recursiva ortogonal), semana 5.
 *
 * Cada nodo cubre una caja del espacio y un tramo contiguo de ranks. Un nodo
 * interno parte su caja por un plano perpendicular a axis y reparte sus ranks
 * entre los dos hijos; una hoja (nprocs == 1) es el dominio de un proceso.
 *
 * El árbol tiene 2P-1 nodos y está REPLICADO en todos los procesos: son unos
 * cientos de bytes, así que cada proceso conoce la caja de todos los demás sin
 * comunicación. Eso es justamente lo que el LET necesita.
 */
typedef struct {
    double min[3], max[3];  /* caja de este nodo */
    int    first, nprocs;   /* tramo de ranks [first, first+nprocs) */
    int    axis;            /* eje de corte; -1 si es hoja */
    double split;           /* coordenada del corte */
    int    left, right;     /* índices de hijos; -1 si es hoja */
} OrbNode;

typedef struct {
    int       rank;
    int       nprocs;
    double    gmin[3], gmax[3];  /* bounding box global cúbico (congelado) */
    uint64_t *splitters;         /* camino Morton: nprocs+1 cortes, el rango de
                                    k es [splitters[k], splitters[k+1]) */
    OrbNode  *orb;               /* camino ORB: 2*nprocs-1 nodos */
    int       orb_count;         /* nodos usados del árbol ORB */
    int       use_orb;           /* 1 = propiedad por caja, 0 = por clave Morton */
} Domain;

/* Reserva el arreglo de splitters y consulta rank/nprocs del comunicador. */
void domain_init(Domain *d, MPI_Comm comm);

/* Libera el arreglo de splitters. */
void domain_free(Domain *d);

/*
 * Calcula el bounding box global de todas las partículas (Allreduce MIN/MAX)
 * y lo convierte en un cubo con margen del 1%, igual que octree_bounds. Todos
 * los procesos obtienen el mismo cubo, condición necesaria para que las claves
 * Morton sean comparables entre procesos.
 *
 * El box se "congela" en d->gmin/gmax y se reutiliza entre reparticiones: si se
 * recodificaran las claves con un box distinto al usado para calcular los
 * splitters, la asignación de dueños sería inconsistente y se perderían
 * partículas. morton_encode acota las coordenadas que se salgan del box.
 */
void domain_global_bounds(Domain *d, const Particle *p, int n_local, MPI_Comm comm);

/*
 * Calcula los nprocs+1 splitters que reparten las n_global partículas en tramos
 * de ~N/P. Colectivo: todos los procesos deben llamarlo y todos obtienen los
 * mismos splitters.
 *
 * Precondición: el arreglo local está ordenado por clave Morton creciente y las
 * claves fueron calculadas con d->gmin/gmax (la búsqueda binaria interna lo
 * asume). Ver morton_sort.
 */
void domain_partition(Domain *d, const Particle *p, int n_local,
                      int64_t n_global, MPI_Comm comm);

/*
 * Igual que domain_partition, pero si use_work != 0 reparte el TRABAJO
 * (Σ p[i].work) en lugar del conteo de partículas.
 *
 * La semana 3 midió que el desbalance de trabajo (1.1149) es bastante mayor que
 * el de partículas (1.0000): domain_partition reparte los conteos casi
 * perfectamente, pero el costo del recorrido no es uniforme y el proceso que
 * recibe la zona densa termina haciendo 11% más trabajo que el promedio. Con
 * pesos se reparte lo que realmente cuesta.
 *
 * El algoritmo no cambia — es la misma bisección simultánea, con la misma
 * cantidad de colectivos — solo cambia la magnitud que se acumula: en vez de
 * "cuántas partículas hay debajo de esta clave", "cuánto trabajo hay debajo de
 * esta clave". Como el arreglo local ya está ordenado por clave, alcanza con un
 * prefijo de work calculado una vez por llamada.
 *
 * El peso está retrasado un paso (es el costo medido en el paso anterior), lo
 * que es válido porque la distribución cambia poco entre pasos consecutivos. En
 * el primer paso work vale 0 para todas: se detecta que el trabajo global es
 * nulo y se cae al reparto por conteo.
 *
 * Efecto esperado y contraintuitivo: el desbalance de PARTÍCULAS empeora (los
 * procesos con zonas densas reciben menos partículas). Es la señal de que
 * funciona, no un problema; la métrica de éxito es el desbalance de force_time.
 */
void domain_partition_ex(Domain *d, const Particle *p, int n_local,
                         int64_t n_global, int use_work, MPI_Comm comm);

/*
 * Particionamiento ORB: bisección recursiva ortogonal del espacio.
 *
 * POR QUÉ EXISTE. La semana 4 midió que el LET no comprime sobre Plummer, y la
 * causa es que un tramo contiguo de claves Morton NO es una región compacta: el
 * AABB del dominio de un proceso llegó a cubrir el 44,6% del espacio, con lo que
 * ningún nodo satisfacía el criterio de exportación y se enviaba todo sin
 * resumir. El problema de fondo es que en Plummer el pico de densidad cae donde
 * se tocan los octantes de nivel 1, así que las claves del core se reparten por
 * todo el espacio de claves y los dominios se interpenetran.
 *
 * ORB corta el espacio por planos: cada proceso queda con una CAJA, los dominios
 * son disjuntos y contiguos, y la caja de un proceso periférico no contiene el
 * core. Recién ahí el criterio del LET puede aceptar resúmenes.
 *
 * ALGORITMO. Recursivo sobre el grupo de procesos: se elige como eje de corte la
 * dimensión más larga de la caja (evita cajas con forma de lámina, que tienen
 * mucha superficie y por lo tanto mucho LET), se bisecta la coordenada hasta que
 * el trabajo quede repartido nL/nprocs, y se recurre en las dos mitades.
 *
 * Se bisecta el TRABAJO (Σ work), no el conteo: el campo work de la semana 4 ya
 * existe y ya está validado, así que ORB hereda el balance por costo sin código
 * nuevo. Con work global nulo (primer paso) cae al conteo.
 *
 * Con P que no es potencia de 2 se parte en nprocs/2 y el resto, con objetivo
 * proporcional. Colectivo.
 */
void domain_partition_orb(Domain *d, const Particle *p, int n_local,
                          int use_work, MPI_Comm comm);

/*
 * Proceso dueño de una posición: desciende el árbol ORB comparando pos[axis]
 * contra el plano de corte. O(log P).
 *
 * Usa el mismo predicado (<) que el conteo del particionamiento. Una discrepancia
 * entre ambos haría que una partícula sobre el plano se contara de un lado y se
 * enviara al otro, y se perderían partículas.
 *
 * No hace test de contención: cualquier posición, aun fuera de la caja global,
 * obtiene un dueño válido. Es lo que hace seguro el caso de una partícula que se
 * escapó del box congelado entre reparticiones.
 */
int domain_owner_pos(const Domain *d, const double pos[3]);

/*
 * Caja del dominio del proceso rank. Con ORB cada proceso conoce la de todos los
 * demás sin comunicación, porque el árbol está replicado.
 */
void domain_box(const Domain *d, int rank, double min[3], double max[3]);

/*
 * Proceso dueño de una clave: el mayor k tal que splitters[k] <= key.
 * Búsqueda binaria sobre los splitters, O(log P).
 */
int domain_owner(const Domain *d, uint64_t key);

#endif /* USE_MPI */

#endif
