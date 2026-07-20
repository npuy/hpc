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
 * Construye el octree de 1000 partículas Plummer y verifica que la masa de la
 * raíz coincida con la suma directa de masas (error relativo < 10⁻¹²). Valida
 * la acumulación bottom-up de octree_compute_mass.
 */
int test_tree_mass(void);

/*
 * Verifica que el centro de masa de la raíz del octree coincida con el CM
 * calculado por suma directa (‖diferencia‖ < 10⁻¹⁰). Valida el promedio
 * ponderado por masa en los nodos internos.
 */
int test_tree_cm(void);

/*
 * Compara la aceleración Barnes-Hut (θ=0.5) contra el método directo O(N²)
 * sobre 1000 partículas Plummer y verifica que el error relativo L2 < 1%.
 * Valida el criterio de apertura y el recorrido del árbol.
 */
int test_bh_force_error(void);

/*
 * Barre θ ∈ {0.8, 0.5, 0.2, 0.0} y verifica que el error respecto al método
 * directo decrezca monótonamente y que θ=0 lo reproduzca (error < 10⁻⁹).
 * Valida la consistencia del esquema de aproximación.
 */
int test_bh_theta_convergence(void);

/*
 * Calcula las fuerzas con 1, 2, 4 y 8 hilos OpenMP y verifica que los tres
 * últimos resultados sean BIT A BIT idénticos al de 1 hilo (memcmp), tanto en
 * compute_forces_bh como en compute_forces_direct_par. Cada iteración escribe
 * solo su propia acc[] y el orden de las sumas no depende del reparto entre
 * hilos, así que cualquier diferencia delata una carrera de datos.
 * Sin OpenMP el test no aplica y se reporta como SKIP.
 */
int test_openmp_determinism(void);

#ifdef USE_MPI

/*
 * Integra 5000 partículas durante 100 pasos con migración y verifica que la
 * suma de n_local sobre todos los procesos sea exactamente N en cada paso.
 * Colectivo: todos los procesos deben invocarlo.
 */
int test_mpi_particle_conservation(void);

/*
 * Verifica que los checksums globales Σ id y Σ id² sean idénticos antes y
 * después de 100 pasos con migración. Detecta pérdidas y duplicados; el
 * segundo momento es necesario porque Σ id sola no distingue un intercambio de
 * dos identidades. Colectivo.
 */
int test_mpi_identity_checksum(void);

/*
 * Verifica que los splitters sean no decrecientes, cubran todo el espacio de
 * claves, que cada partícula local caiga dentro del tramo de su proceso, y que
 * el desbalance max/avg del reparto sea menor a 1.15. Colectivo.
 */
int test_mpi_partition_validity(void);

/*
 * Compara las aceleraciones del camino distribuido completo (particionar,
 * migrar, replicar, calcular el tramo propio) contra el resultado del cálculo
 * secuencial sobre el mismo estado inicial, emparejando las partículas por id
 * porque la migración las permuta. Criterio: error relativo máximo < 1e-12.
 * Colectivo.
 */
int test_mpi_forces_vs_sequential(void);

#endif /* USE_MPI */

/*
 * Ejecuta la suite completa e imprime un resumen. Son 9 tests en la versión
 * secuencial y 13 en la versión MPI (los 4 distribuidos requieren mpirun).
 * Retorna 0 si todos pasan, 1 si alguno falla (compatible con exit code).
 */
int run_all_tests(void);

#endif
