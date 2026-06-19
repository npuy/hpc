#ifndef VALIDATION_H
#define VALIDATION_H

/*
 * validation.h — Suite de tests de validación física y numérica.
 *
 * Verifica propiedades conservadas del sistema (energía, momento lineal),
 * correctitud del integrador (órbita circular analítica) y del ordenamiento
 * espacial (Morton). Cada test retorna 1 si pasa, 0 si falla.
 * Se invoca desde main con la flag --validate.
 */

/*
 * Integra un sistema de 2 cuerpos en órbita circular durante un período
 * completo (10.000 pasos) y verifica que las posiciones finales coincidan
 * con las iniciales con error relativo < 1%. Valida conjuntamente el
 * integrador Leapfrog y el cálculo de fuerzas.
 */
int test_two_body_orbit(void);

/*
 * Integra 100 partículas Plummer durante 100 pasos y verifica que el drift
 * relativo de energía total |ΔE/E₀| < 10⁻³. Confirma que el integrador
 * simpléctico no introduce drift secular en la energía.
 */
int test_energy_conservation(void);

/*
 * Integra 50 partículas Plummer durante 50 pasos y verifica que el cambio
 * en momento lineal total |ΔP| < 10⁻¹². Valida la simetría de Newton
 * (acción = -reacción) en compute_forces_direct.
 */
int test_momentum_conservation(void);

/*
 * Codifica las 8 esquinas del cubo unitario [0,1]³ con Morton y verifica
 * que morton_sort produzca un orden estrictamente creciente con 8 claves
 * distintas. Valida la correctitud del bit interleaving.
 */
int test_morton_ordering(void);

/*
 * Ejecuta los 4 tests secuencialmente e imprime un resumen de resultados.
 * Retorna 0 si todos pasan, 1 si alguno falla (compatible con exit code).
 */
int run_all_tests(void);

#endif
