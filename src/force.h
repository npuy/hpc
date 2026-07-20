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
 * Variante paralelizable con OpenMP del método directo. NO usa la tercera ley
 * de Newton: cada partícula i acumula sobre una copia local recorriendo todos
 * los j != i, y escribe solo en p[i].acc. Cuesta el doble de flops que
 * compute_forces_direct pero no tiene carreras de datos y es bit a bit
 * reproducible con cualquier número de hilos.
 *
 * compute_forces_direct escribe en p[j].acc dentro del bucle interno y por eso
 * NO puede paralelizarse sobre i; se conserva como oráculo secuencial exacto.
 * Los resultados de ambas difieren solo en el orden de las sumas.
 */
void compute_forces_direct_par(Particle *p, int n, double softening);

/*
 * Energía cinética total del sistema: K = Σᵢ ½ mᵢ |vᵢ|².
 * Necesaria junto con potential_energy para monitorear conservación de energía.
 */
double kinetic_energy(const Particle *p, int n);

/*
 * Contribución de las partículas [i0, i1) a la energía cinética total.
 * Sumando el resultado de todos los procesos MPI (Allreduce) se obtiene la
 * energía cinética global. kinetic_energy es el caso particular (0, n).
 */
double kinetic_energy_range(const Particle *p, int n, int i0, int i1);

/*
 * Energía potencial gravitatoria total: U = -Σᵢ<ⱼ mᵢmⱼ / √(|rᵢⱼ|² + ε²).
 * Usa el mismo softening que compute_forces_direct para consistencia del
 * hamiltoniano: E = K + U debe ser integral de movimiento del sistema integrado.
 */
double potential_energy(const Particle *p, int n, double softening);

#endif
