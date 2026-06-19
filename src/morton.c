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

    uint64_t ix = (uint64_t)(((x - min[0]) / dx) * range);
    uint64_t iy = (uint64_t)(((y - min[1]) / dy) * range);
    uint64_t iz = (uint64_t)(((z - min[2]) / dz) * range);

    return (expand_bits(ix) << 2) | (expand_bits(iy) << 1) | expand_bits(iz);
}

static int morton_compare(const void *a, const void *b) {
    const Particle *pa = (const Particle *)a;
    const Particle *pb = (const Particle *)b;
    if (pa->morton < pb->morton) return -1;
    if (pa->morton > pb->morton) return 1;
    return 0;
}

void morton_sort(Particle *p, int n, double min[3], double max[3]) {
    for (int i = 0; i < n; i++) {
        p[i].morton = morton_encode(p[i].pos[0], p[i].pos[1], p[i].pos[2],
                                    min, max, 21);
    }
    qsort(p, n, sizeof(Particle), morton_compare);
}
