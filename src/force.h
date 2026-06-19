#ifndef FORCE_H
#define FORCE_H

/*
 * force.h — Cálculo de interacciones gravitatorias y observables energéticos.
 *
 * Implementa el método directo (all-pairs) O(N²) con simetría de Newton:
 * cada par (i,j) se evalúa una sola vez y la fuerza se acumula en ambas
 * partículas con signo opuesto. Esto reduce las evaluaciones a N(N-1)/2
 * y garantiza conservación exacta del momento lineal.
 *
 * El softening ε regulariza la fuerza en distancias cortas sumando ε²
 * a |r|² antes de invertir, evitando aceleraciones divergentes.
 */

#include "particle.h"

/*
 * Calcula la aceleración gravitatoria de cada partícula por interacción directa
 * con todas las demás. Resetea acc[] a cero antes de acumular. Usa la tercera
 * ley de Newton para reducir evaluaciones. Complejidad: O(N²) tiempo, O(1) espacio.
 */
void compute_forces_direct(Particle *p, int n, double softening);

/*
 * Energía cinética total del sistema: K = Σᵢ ½ mᵢ |vᵢ|².
 * Necesaria junto con potential_energy para monitorear conservación de energía.
 */
double kinetic_energy(const Particle *p, int n);

/*
 * Energía potencial gravitatoria total: U = -Σᵢ<ⱼ mᵢmⱼ / √(|rᵢⱼ|² + ε²).
 * Usa el mismo softening que compute_forces_direct para consistencia del
 * hamiltoniano: E = K + U debe ser integral de movimiento del sistema integrado.
 */
double potential_energy(const Particle *p, int n, double softening);

#endif
