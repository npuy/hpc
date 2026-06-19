#ifndef VEC3_H
#define VEC3_H

/*
 * vec3.h — Aritmética vectorial 3D mediante macros.
 *
 * Todas las macros operan sobre arreglos double[3].
 * Las macros-sentencia usan el patrón do{...}while(0) para ser seguras en
 * contextos de control de flujo. Las macros-expresión (DOT, NORM) retornan
 * un double y pueden usarse en asignaciones o condicionales.
 *
 * ADVERTENCIA: los argumentos se evalúan múltiples veces; no pasar
 * expresiones con efectos secundarios (ej: a++).
 */

#include <math.h>

/* Resetea las 3 componentes de a a cero */
#define VEC3_ZERO(a) do { (a)[0]=0; (a)[1]=0; (a)[2]=0; } while(0)
/* Copia componente a componente: r ← a */
#define VEC3_COPY(r, a) do { (r)[0]=(a)[0]; (r)[1]=(a)[1]; (r)[2]=(a)[2]; } while(0)
/* Suma vectorial: r ← a + b */
#define VEC3_ADD(r, a, b) do { (r)[0]=(a)[0]+(b)[0]; (r)[1]=(a)[1]+(b)[1]; (r)[2]=(a)[2]+(b)[2]; } while(0)
/* Resta vectorial: r ← a - b */
#define VEC3_SUB(r, a, b) do { (r)[0]=(a)[0]-(b)[0]; (r)[1]=(a)[1]-(b)[1]; (r)[2]=(a)[2]-(b)[2]; } while(0)
/* Escalado: r ← a * s, donde s es un escalar double */
#define VEC3_SCALE(r, a, s) do { (r)[0]=(a)[0]*(s); (r)[1]=(a)[1]*(s); (r)[2]=(a)[2]*(s); } while(0)
/* Multiply-add: r ← a + b * s (usado en integración temporal: pos += vel * dt) */
#define VEC3_MADD(r, a, b, s) do { (r)[0]=(a)[0]+(b)[0]*(s); (r)[1]=(a)[1]+(b)[1]*(s); (r)[2]=(a)[2]+(b)[2]*(s); } while(0)
/* Producto punto: retorna a·b = Σ aᵢbᵢ */
#define VEC3_DOT(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
/* Norma euclidiana: retorna ‖a‖ = √(a·a) */
#define VEC3_NORM(a) sqrt(VEC3_DOT(a, a))

#endif
