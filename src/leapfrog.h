#ifndef LEAPFROG_H
#define LEAPFROG_H

/*
 * leapfrog.h — Integrador simpléctico Leapfrog (variante KDK).
 *
 * Provee los dos operadores elementales del esquema Kick-Drift-Kick:
 *   1. Kick  (½dt): v += a * dt/2
 *   2. Drift (dt):  x += v * dt
 *
 * Un paso completo se construye como: Kick → Drift → [recalcular fuerzas] → Kick.
 * El esquema es de segundo orden, simpléctico (sin drift secular en energía),
 * time-reversible, y conserva el momento lineal exactamente.
 */

#include "particle.h"

/*
 * Aplica un medio-paso de velocidad: vᵢ += aᵢ * dt_half para todas las
 * partículas. Requiere que acc[] haya sido previamente calculado por
 * compute_forces_direct. El argumento dt_half es dt/2, no dt completo.
 */
void leapfrog_kick(Particle *p, int n, double dt_half);

/*
 * Aplica un paso completo de posición: xᵢ += vᵢ * dt para todas las
 * partículas. Las velocidades deben estar en el punto medio temporal
 * (después de un kick) para que el esquema sea de segundo orden.
 */
void leapfrog_drift(Particle *p, int n, double dt);

#endif
