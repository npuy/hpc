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

typedef struct {
    int       rank;
    int       nprocs;
    double    gmin[3], gmax[3];  /* bounding box global cúbico (congelado) */
    uint64_t *splitters;         /* nprocs+1 cortes: el rango de k es
                                    [splitters[k], splitters[k+1]) */
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
 * Proceso dueño de una clave: el mayor k tal que splitters[k] <= key.
 * Búsqueda binaria sobre los splitters, O(log P).
 */
int domain_owner(const Domain *d, uint64_t key);

#endif /* USE_MPI */

#endif
