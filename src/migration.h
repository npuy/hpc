#ifndef MIGRATION_H
#define MIGRATION_H

/*
 * migration.h — Reubicación de partículas entre procesos (etapa 7).
 *
 * Después del drift las partículas se movieron y algunas quedaron fuera del
 * tramo de claves Morton de su proceso. migrate_particles recalcula las claves,
 * determina el dueño de cada partícula y hace el intercambio con MPI_Alltoallv.
 *
 * Se usa Alltoallv y no un patrón Send/Recv armado a mano por dos razones: es
 * colectivo, así que no puede deadlockear por un orden de envíos mal elegido, y
 * el runtime elige el algoritmo interno (bruck para mensajes chicos, pairwise
 * para grandes) según el tamaño.
 *
 * Este archivo contiene además replicate_particles, la réplica global que la
 * semana 3 usa en lugar del LET. Está aislada en una sola función a propósito:
 * la semana 4 la sustituye por el intercambio de árboles localmente esenciales
 * tocando únicamente este punto.
 */

#ifdef USE_MPI

#include <mpi.h>
#include "particle.h"
#include "domain.h"

/*
 * Reubica las partículas cuya clave salió del tramo local. Recalcula las claves
 * Morton con el bounding box global congelado (d->gmin/gmax), reparte por
 * destino y hace el intercambio. Al volver, el arreglo local queda ordenado por
 * clave creciente, que es la invariante que asume domain_partition.
 *
 * *p puede ser reasignado (realloc) si llegan más partículas de las que caben;
 * *capacity se actualiza. Retorna el nuevo n_local, que puede ser 0.
 *
 * Colectivo: todos los procesos del comunicador deben llamarlo.
 */
int migrate_particles(Particle **p, int n_local, int *capacity,
                      const Domain *d, MPI_Datatype ptype, MPI_Comm comm);

/*
 * Reúne en all[] las partículas de todos los procesos (MPI_Allgatherv) y
 * devuelve en *offset el índice donde empiezan las locales dentro de all[].
 *
 * Como cada proceso mantiene sus partículas ordenadas por clave y los tramos de
 * claves son disjuntos y crecientes con el rank, la concatenación que produce
 * Allgatherv YA está globalmente ordenada por Morton: no hace falta ningún sort
 * posterior, y las partículas locales ocupan el rango contiguo
 * [*offset, *offset + n_local) — que es justo lo que compute_forces_bh_range
 * necesita.
 *
 * PROVISIONAL (semana 3): cuesta O(N) en memoria y comunicación por proceso, y
 * es el techo de escalabilidad de esta etapa. La semana 4 lo reemplaza por LET.
 *
 * all[] debe tener capacidad para n_global partículas. Colectivo.
 */
void replicate_particles(const Particle *local, int n_local,
                         Particle *all, int n_global, int *offset,
                         MPI_Datatype ptype, MPI_Comm comm);

/*
 * Checksums globales de identidad: *sum_id = Σ id y *sum_id2 = Σ id².
 * Si ambos coinciden con sus valores iniciales, no se perdió ni se duplicó
 * ninguna partícula. El segundo momento es necesario: Σ id sola no detecta que
 * una partícula haya sido reemplazada por otra de id complementario.
 * Colectivo.
 */
void particle_checksum(const Particle *p, int n_local,
                       int64_t *sum_id, int64_t *sum_id2, MPI_Comm comm);

#endif /* USE_MPI */

#endif
