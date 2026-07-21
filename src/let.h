#ifndef LET_H
#define LET_H

/*
 * let.h — Árbol localmente esencial (Locally Essential Tree), etapa 8.
 *
 * Reemplaza la réplica global de la semana 3. En vez de que cada proceso reciba
 * las N partículas y construya el árbol completo, cada proceso recibe solo el
 * subconjunto del árbol global que necesita para calcular las fuerzas de SUS
 * partículas con la misma precisión que tendría con el árbol entero. Todo lo
 * que está suficientemente lejos llega resumido en unos pocos nodos agregados.
 *
 * Por qué importa: la semana 3 midió que el Allgatherv era apenas el 1.3% del
 * tiempo, así que el LET NO se justifica por ahorro de red. Lo caro era la
 * construcción del árbol replicada — tree_time crecía 4.8× de P=1 a P=8 porque
 * cada proceso construía sobre las N partículas. El LET elimina ese trabajo
 * redundante: cada proceso construye sobre sus N/P locales más los fantasmas.
 *
 * CRITERIO DE EXPORTACIÓN (MAC conservador). Para el nodo c de lado s y centro
 * de masa cm, y la caja destino B_r:
 *
 *     d_min = distancia mínima de cm a B_r   (0 si cm cae dentro de B_r)
 *     aceptar el resumen si   s² < theta² · (d_min² + eps²)
 *
 * Es correcto porque el receptor, al calcular la fuerza sobre su partícula
 * i ∈ B_r, aplicaría s² < theta²·(d_i² + eps²) con d_i ≥ d_min. Usar la
 * distancia MÍNIMA sobre toda la caja garantiza que nunca exportamos un resumen
 * donde el receptor quería detalle. El caso inverso (abrir un nodo que el
 * receptor habría aceptado) manda datos de más: cuesta ancho de banda, no
 * precisión.
 *
 * REPRESENTACIÓN. Lo exportado viaja como Particle: los nodos aceptados como
 * pseudo-partículas (pos = cm, mass = masa del nodo, id = PARTICLE_GHOST_ID) y
 * las hojas como sus partículas reales. Así se reutilizan el tipo MPI derivado
 * y MPI_Alltoallv tal cual, sin un tipo nuevo ni empaquetado a mano.
 *
 * SE PIERDE LA EXACTITUD BIT A BIT, y es correcto que así sea. Un fantasma
 * llega como punto de masa en cm; al reinsertarlo en el árbol puede quedar
 * agrupado con otros bajo un ancestro común que a su vez satisface el MAC, y se
 * aproxima dos veces. El error sigue acotado por theta pero ya no es cero, así
 * que el test 13 de la semana 3 (error exactamente 0.0) no se puede mantener
 * como criterio. Lo reemplaza el test 14: con theta = 0 nada se aproxima, el
 * LET degenera en réplica completa y la igualdad vuelve a ser bit a bit — es el
 * test que detecta cualquier interacción faltante.
 */

#ifdef USE_MPI

#include <mpi.h>
#include "particle.h"
#include "octree.h"
#include "domain.h"

/*
 * Publica el AABB de las partículas locales y recibe el de todos los procesos
 * en boxes[6*P] con el formato (min[3], max[3]) por proceso.
 *
 * Hay que llamarlo CADA PASO, después de migrar y antes de exportar: las
 * partículas se movieron y una caja obsoleta puede dejar afuera interacciones
 * necesarias. El costo es un Allgather de 6 doubles por proceso, despreciable
 * frente al Alltoallv del propio LET.
 *
 * Un proceso sin partículas publica una caja degenerada en el centro del
 * dominio: no necesita fantasmas (no calcula ninguna fuerza), y lo que le
 * manden los demás es inofensivo.
 */
void let_gather_domain_boxes(const Particle *local, int n_local,
                             const Domain *d, double *boxes, MPI_Comm comm);

/*
 * Selecciona, intercambia e importa el árbol localmente esencial.
 *
 * Recorre local_tree una vez por destino aplicando el MAC conservador, arma un
 * paquete por proceso y hace el intercambio con MPI_Alltoallv. Los fantasmas
 * recibidos quedan APPENDEADOS en *p después de las n_local locales, en el
 * rango [n_local, n_local + n_ghost). *p puede reasignarse (realloc) y
 * *capacity se actualiza. Retorna n_ghost.
 *
 * Precondiciones: local_tree fue construido con octree_build_box sobre la caja
 * de d (todos los procesos, la misma caja) y ya tiene masas y centros de masa
 * calculados — el criterio de exportación los necesita.
 *
 * Los fantasmas son transitorios por diseño: el migrate_particles del paso
 * siguiente reemplaza el arreglo local y los descarta sin que haya que
 * limpiarlos explícitamente. Quien recorra el arreglo debe usar el rango
 * [0, n_local) para todo lo que sea checksums, energías o migración: un
 * fantasma tiene id = -1 y masa de un nodo agregado, así que contarlo rompe el
 * checksum e infla la energía.
 *
 * Colectivo: todos los procesos del comunicador deben llamarlo.
 */
int let_exchange(Particle **p, int n_local, int *capacity,
                 const Octree *local_tree, const Domain *d,
                 double theta, double softening,
                 MPI_Datatype ptype, MPI_Comm comm);

/*
 * Verificación barata e independiente de que la selección no descartó masa:
 * Σ mass sobre locales + fantasmas de un proceso debe dar la masa total del
 * sistema, porque el LET resume pero no tira nada. Retorna el error relativo
 * máximo entre procesos.
 *
 * Es la primera comprobación a mirar cuando las fuerzas salen mal: detecta casi
 * cualquier bug de selección (podas de más, nodos duplicados, hojas vacías mal
 * manejadas) sin necesidad de comparar contra una referencia secuencial.
 */
double let_mass_error(const Particle *p, int n_local, int n_ghost,
                      MPI_Comm comm);

#endif /* USE_MPI */

#endif
