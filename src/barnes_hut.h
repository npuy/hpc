#ifndef BARNES_HUT_H
#define BARNES_HUT_H

/*
 * barnes_hut.h — Cálculo de fuerzas O(N log N) por recorrido del octree.
 *
 * Para cada partícula se recorre el árbol desde la raíz aplicando el criterio
 * de apertura (MAC, Multipole Acceptance Criterion): si una celda de lado s
 * está a distancia d y cumple s/d < θ, todo su contenido se aproxima por su
 * masa y centro de masa como un único cuerpo; en caso contrario se desciende a
 * los 8 hijos. Con θ=0 no se aproxima nada y el resultado coincide con el
 * método directo O(N²) (ver force.h); θ≈0.5 es el estándar.
 */

#include "particle.h"
#include "octree.h"

/*
 * Calcula la aceleración de cada partícula recorriendo el octree t con ángulo
 * de apertura theta y softening. Resetea acc[] antes de acumular. Usa la misma
 * fórmula gravitatoria con softening que compute_forces_direct, de modo que
 * theta=0 reproduce el método directo salvo orden de sumas. El árbol debe
 * tener masas/centros de masa ya calculados (octree_compute_mass).
 */
void compute_forces_bh(const Octree *t, Particle *p, int n,
                       double theta, double softening);

/*
 * Energía potencial gravitatoria total aproximada por Barnes-Hut:
 * U = ½ Σᵢ mᵢ φᵢ, con φᵢ el potencial en la partícula i evaluado por el mismo
 * recorrido del árbol. Permite monitorear conservación de energía en O(N log N)
 * en vez del O(N²) de potential_energy. Con theta=0 coincide con el directo.
 */
double potential_energy_bh(const Octree *t, const Particle *p, int n,
                           double theta, double softening);

#endif
