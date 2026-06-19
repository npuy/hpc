# Módulo: leapfrog

## Descripción general

Implementa el integrador simpléctico Leapfrog en su variante KDK (Kick-Drift-Kick). Este esquema es estándar en simulaciones gravitatorias de N cuerpos porque conserva la estructura simpléctica del espacio de fases, lo que garantiza estabilidad a largo plazo en la energía total (sin drift secular) y conservación exacta del momento lineal.

## Esquema de integración KDK

Un paso temporal completo se compone de tres sub-pasos:

```
1. Kick   (½ dt):   v(t)      → v(t + ½dt)   usando a(t)
2. Drift  (dt):     x(t)      → x(t + dt)     usando v(t + ½dt)
3. [Recalcular fuerzas en x(t + dt)]
4. Kick   (½ dt):   v(t + ½dt) → v(t + dt)    usando a(t + dt)
```

La variante KDK tiene la ventaja de que posiciones y velocidades quedan sincronizadas al mismo tiempo `t`, lo cual facilita el cálculo de observables (energía, momento).

## Interfaz pública

### `leapfrog_kick(p, n, dt_half)`

Aplica un medio-paso de velocidad a todas las `n` partículas: `vᵢ += aᵢ * dt/2`. Requiere que las aceleraciones `acc[]` hayan sido previamente calculadas por `compute_forces_direct` (o su equivalente).

**Parámetros**:
- `p`: arreglo de partículas (lectura/escritura de `vel[]`, lectura de `acc[]`)
- `n`: cantidad de partículas
- `dt_half`: **medio** paso temporal (`dt/2`), no el paso completo

### `leapfrog_drift(p, n, dt)`

Aplica un paso completo de posición a todas las `n` partículas: `xᵢ += vᵢ * dt`.

**Parámetros**:
- `p`: arreglo de partículas (lectura/escritura de `pos[]`, lectura de `vel[]`)
- `n`: cantidad de partículas
- `dt`: paso temporal completo

## Dependencias

- `particle.h` — tipo `Particle`.
- `vec3.h` — incluido pero no usado directamente (las operaciones se hacen componente a componente para claridad).

## Propiedades del integrador

| Propiedad | Valor |
|---|---|
| Orden de convergencia | 2 (error global O(dt²)) |
| Simpléctico | Sí |
| Conservación de energía | Oscila sin drift secular |
| Conservación de momento | Exacta (a precisión de máquina) |
| Time-reversible | Sí |

## Protocolo de uso

```c
compute_forces_direct(p, n, eps);   // a(t=0)
for (t = 0; t < t_end; t += dt) {
    leapfrog_kick(p, n, dt / 2.0);   // v += a * dt/2
    leapfrog_drift(p, n, dt);         // x += v * dt
    compute_forces_direct(p, n, eps); // a(t+dt)
    leapfrog_kick(p, n, dt / 2.0);   // v += a * dt/2
}
```
