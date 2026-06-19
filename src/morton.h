#ifndef MORTON_H
#define MORTON_H

/*
 * morton.h — Codificación Morton (Z-order curve) para ordenamiento espacial.
 *
 * Mapea coordenadas 3D continuas a claves enteras de 63 bits intercalando
 * los bits de cada eje (x en bits 2,5,8,...; y en 1,4,7,...; z en 0,3,6,...).
 * Partículas cercanas en el espacio tienden a tener claves cercanas, lo que
 * al ordenar por clave agrupa vecinos espaciales en memoria contigua.
 * Esto mejora la coherencia de caché en el cálculo de fuerzas y es la base
 * para la futura construcción de octrees (Barnes-Hut).
 */

#include "particle.h"
#include <stdint.h>

/*
 * Codifica una posición (x,y,z) continua en una clave Morton entera.
 * Normaliza las coordenadas al rango [0, 2^bits - 1] usando el bounding
 * box [min, max], discretiza a enteros, e intercala los bits de los 3 ejes.
 * Con bits=21 produce claves de 63 bits (~2M celdas por dimensión).
 */
uint64_t morton_encode(double x, double y, double z,
                       double min[3], double max[3], int bits);

/*
 * Calcula claves Morton con 21 bits/eje para todas las n partículas según
 * su posición dentro del bounding box [min, max], y luego ordena el arreglo
 * in-place por clave creciente (qsort). Modifica tanto el campo morton
 * como el orden de los elementos en el arreglo.
 */
void morton_sort(Particle *p, int n, double min[3], double max[3]);

#endif
