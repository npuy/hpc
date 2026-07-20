#include "morton.h"
#include <stdlib.h>

static uint64_t expand_bits(uint64_t v) {
    v &= 0x1fffff;
    v = (v | v << 32) & 0x1f00000000ffffULL;
    v = (v | v << 16) & 0x1f0000ff0000ffULL;
    v = (v | v << 8)  & 0x100f00f00f00f00fULL;
    v = (v | v << 4)  & 0x10c30c30c30c30c3ULL;
    v = (v | v << 2)  & 0x1249249249249249ULL;
    return v;
}

uint64_t morton_encode(double x, double y, double z,
                       double min[3], double max[3], int bits) {
    double range = (1ULL << bits) - 1;
    double dx = max[0] - min[0];
    double dy = max[1] - min[1];
    double dz = max[2] - min[2];
    if (dx < 1e-15) dx = 1e-15;
    if (dy < 1e-15) dy = 1e-15;
    if (dz < 1e-15) dz = 1e-15;

    /*
     * Acotar antes del cast. Con el bounding box congelado entre reparticiones
     * (semana 3) una particula puede salirse del box: sin clamp, un valor > range
     * se envuelve por el &= 0x1fffff de expand_bits (el extremo derecho recibiria
     * una clave del extremo izquierdo y migraria al proceso equivocado), y un
     * valor negativo convertido a uint64_t es comportamiento indefinido.
     */
    double fx = ((x - min[0]) / dx) * range;
    double fy = ((y - min[1]) / dy) * range;
    double fz = ((z - min[2]) / dz) * range;
    if (fx < 0.0) fx = 0.0; else if (fx > range) fx = range;
    if (fy < 0.0) fy = 0.0; else if (fy > range) fy = range;
    if (fz < 0.0) fz = 0.0; else if (fz > range) fz = range;

    uint64_t ix = (uint64_t)fx;
    uint64_t iy = (uint64_t)fy;
    uint64_t iz = (uint64_t)fz;

    return (expand_bits(ix) << 2) | (expand_bits(iy) << 1) | expand_bits(iz);
}

static int morton_compare(const void *a, const void *b) {
    const Particle *pa = (const Particle *)a;
    const Particle *pb = (const Particle *)b;
    if (pa->morton < pb->morton) return -1;
    if (pa->morton > pb->morton) return 1;
    return 0;
}

void morton_keys(Particle *p, int n, double min[3], double max[3]) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        p[i].morton = morton_encode(p[i].pos[0], p[i].pos[1], p[i].pos[2],
                                    min, max, 21);
    }
}

void morton_sort(Particle *p, int n, double min[3], double max[3]) {
    morton_keys(p, n, min, max);
    /* El qsort queda serial: paralelizarlo requiere un merge sort propio y es
       una fraccion chica del paso comparada con el calculo de fuerzas. Junto
       con octree_build forma la parte serial que acota el speedup (Amdahl). */
    qsort(p, n, sizeof(Particle), morton_compare);
}
