#include "validation.h"
#include "particle.h"
#include "force.h"
#include "leapfrog.h"
#include "morton.h"
#include "octree.h"
#include "barnes_hut.h"
#include "vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int test_two_body_orbit(void) {
    printf("=== Test 1: Orbita circular de dos cuerpos ===\n");

    Particle p[2];
    int n;
    init_two_body(p, &n);

    double pos0_x = p[0].pos[0];
    double pos0_y = p[0].pos[1];
    double pos1_x = p[1].pos[0];
    double pos1_y = p[1].pos[1];

    double softening = 0.001;
    double sep = p[1].pos[0] - p[0].pos[0];
    /* Period = 2*pi*r^(3/2) / sqrt(G*M) with r=sep, G=1, M=1 */
    double T = 2.0 * M_PI * pow(sep, 1.5);
    double dt = T / 10000.0;
    int steps = (int)(T / dt);

    compute_forces_direct(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double err0 = sqrt((p[0].pos[0]-pos0_x)*(p[0].pos[0]-pos0_x) +
                        (p[0].pos[1]-pos0_y)*(p[0].pos[1]-pos0_y));
    double err1 = sqrt((p[1].pos[0]-pos1_x)*(p[1].pos[0]-pos1_x) +
                        (p[1].pos[1]-pos1_y)*(p[1].pos[1]-pos1_y));
    double max_err = fmax(err0, err1);
    double ref = sep;

    printf("  Periodo orbital T = %.4f, dt = %.6f, pasos = %d\n", T, dt, steps);
    printf("  Error pos cuerpo 0: %.6e\n", err0);
    printf("  Error pos cuerpo 1: %.6e\n", err1);
    printf("  Error relativo max: %.4f%%\n", (max_err / ref) * 100.0);

    int pass = (max_err / ref) < 0.01;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    return pass;
}

int test_energy_conservation(void) {
    printf("=== Test 2: Conservacion de energia (N=100, 100 pasos) ===\n");

    int n = 100;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double softening = 0.01;
    double dt = 0.001;
    int steps = 100;

    compute_forces_direct(p, n, softening);
    double E0 = kinetic_energy(p, n) + potential_energy(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double Ef = kinetic_energy(p, n) + potential_energy(p, n, softening);
    double drift = fabs(Ef - E0) / fabs(E0);

    printf("  E0 = %.10e\n", E0);
    printf("  Ef = %.10e\n", Ef);
    printf("  Drift relativo: %.6e\n", drift);

    int pass = drift < 1e-3;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(p);
    return pass;
}

int test_momentum_conservation(void) {
    printf("=== Test 3: Conservacion de momento lineal ===\n");

    int n = 50;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 123);

    double softening = 0.01;
    double dt = 0.001;
    int steps = 50;

    double p0[3] = {0, 0, 0};
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            p0[k] += p[i].mass * p[i].vel[k];

    compute_forces_direct(p, n, softening);

    for (int s = 0; s < steps; s++) {
        leapfrog_kick(p, n, dt / 2.0);
        leapfrog_drift(p, n, dt);
        compute_forces_direct(p, n, softening);
        leapfrog_kick(p, n, dt / 2.0);
    }

    double pf[3] = {0, 0, 0};
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            pf[k] += p[i].mass * p[i].vel[k];

    double dp = sqrt((pf[0]-p0[0])*(pf[0]-p0[0]) +
                      (pf[1]-p0[1])*(pf[1]-p0[1]) +
                      (pf[2]-p0[2])*(pf[2]-p0[2]));

    printf("  Momento inicial: (%.6e, %.6e, %.6e)\n", p0[0], p0[1], p0[2]);
    printf("  Momento final:   (%.6e, %.6e, %.6e)\n", pf[0], pf[1], pf[2]);
    printf("  |dp|: %.6e\n", dp);

    int pass = dp < 1e-12;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(p);
    return pass;
}

int test_morton_ordering(void) {
    printf("=== Test 4: Morton ordering (8 esquinas del cubo) ===\n");

    Particle p[8];
    memset(p, 0, sizeof(p));
    double corners[8][3] = {
        {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
        {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}
    };
    for (int i = 0; i < 8; i++) {
        p[i].pos[0] = corners[i][0];
        p[i].pos[1] = corners[i][1];
        p[i].pos[2] = corners[i][2];
        p[i].mass = 1.0;
    }

    double mn[3] = {0, 0, 0};
    double mx[3] = {1, 1, 1};

    for (int i = 0; i < 8; i++) {
        p[i].morton = morton_encode(p[i].pos[0], p[i].pos[1], p[i].pos[2],
                                    mn, mx, 21);
    }

    printf("  Claves Morton antes de sort:\n");
    for (int i = 0; i < 8; i++)
        printf("    (%.0f,%.0f,%.0f) -> 0x%016llx\n",
               p[i].pos[0], p[i].pos[1], p[i].pos[2],
               (unsigned long long)p[i].morton);

    morton_sort(p, 8, mn, mx);

    printf("  Claves Morton despues de sort:\n");
    for (int i = 0; i < 8; i++)
        printf("    (%.0f,%.0f,%.0f) -> 0x%016llx\n",
               p[i].pos[0], p[i].pos[1], p[i].pos[2],
               (unsigned long long)p[i].morton);

    int pass = 1;
    for (int i = 1; i < 8; i++) {
        if (p[i].morton < p[i-1].morton) {
            pass = 0;
            break;
        }
    }

    /* Verify all 8 distinct keys */
    for (int i = 0; i < 8 && pass; i++)
        for (int j = i+1; j < 8 && pass; j++)
            if (p[i].morton == p[j].morton)
                pass = 0;

    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");
    return pass;
}

/* Error relativo L2 global entre dos campos de aceleracion. */
static double accel_l2_error(const Particle *ref, const Particle *bh, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        double d[3];
        VEC3_SUB(d, bh[i].acc, ref[i].acc);
        num += VEC3_DOT(d, d);
        den += VEC3_DOT(ref[i].acc, ref[i].acc);
    }
    return sqrt(num / den);
}

int test_tree_mass(void) {
    printf("=== Test 5: Conservacion de masa del arbol ===\n");

    int n = 1000;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double total = 0.0;
    for (int i = 0; i < n; i++) total += p[i].mass;

    Octree *t = octree_build(p, n);
    octree_compute_mass(t, p);
    double root_mass = t->nodes[t->root].mass;

    double err = fabs(root_mass - total) / total;
    printf("  Masa total (suma directa): %.15e\n", total);
    printf("  Masa de la raiz:           %.15e\n", root_mass);
    printf("  Error relativo: %.6e\n", err);

    int pass = err < 1e-12;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(p);
    return pass;
}

int test_tree_cm(void) {
    printf("=== Test 6: Centro de masa del arbol ===\n");

    int n = 1000;
    Particle *p = malloc(n * sizeof(Particle));
    init_plummer(p, n, 42);

    double total = 0.0, cm[3] = {0, 0, 0};
    for (int i = 0; i < n; i++) {
        total += p[i].mass;
        for (int k = 0; k < 3; k++) cm[k] += p[i].mass * p[i].pos[k];
    }
    for (int k = 0; k < 3; k++) cm[k] /= total;

    Octree *t = octree_build(p, n);
    octree_compute_mass(t, p);
    double *rcm = t->nodes[t->root].cm;

    double d[3];
    VEC3_SUB(d, rcm, cm);
    double err = VEC3_NORM(d);
    printf("  CM directo: (%.6e, %.6e, %.6e)\n", cm[0], cm[1], cm[2]);
    printf("  CM raiz:    (%.6e, %.6e, %.6e)\n", rcm[0], rcm[1], rcm[2]);
    printf("  |diferencia|: %.6e\n", err);

    int pass = err < 1e-10;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(p);
    return pass;
}

int test_bh_force_error(void) {
    printf("=== Test 7: Error de fuerza Barnes-Hut vs directo (theta=0.5) ===\n");

    int n = 1000;
    double softening = 0.01, theta = 0.5;
    Particle *ref = malloc(n * sizeof(Particle));
    init_plummer(ref, n, 42);

    Particle *bh = malloc(n * sizeof(Particle));
    memcpy(bh, ref, n * sizeof(Particle));

    compute_forces_direct(ref, n, softening);

    Octree *t = octree_build(bh, n);
    octree_compute_mass(t, bh);
    compute_forces_bh(t, bh, n, theta, softening);

    double err = accel_l2_error(ref, bh, n);
    printf("  Error relativo L2 de la aceleracion: %.6e\n", err);

    int pass = err < 1e-2;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    octree_free(t);
    free(bh);
    free(ref);
    return pass;
}

int test_bh_theta_convergence(void) {
    printf("=== Test 8: Convergencia en theta ===\n");

    int n = 1000;
    double softening = 0.01;
    double thetas[4] = {0.8, 0.5, 0.2, 0.0};

    Particle *ref = malloc(n * sizeof(Particle));
    init_plummer(ref, n, 42);
    compute_forces_direct(ref, n, softening);

    double errs[4];
    for (int idx = 0; idx < 4; idx++) {
        Particle *bh = malloc(n * sizeof(Particle));
        memcpy(bh, ref, n * sizeof(Particle));
        Octree *t = octree_build(bh, n);
        octree_compute_mass(t, bh);
        compute_forces_bh(t, bh, n, thetas[idx], softening);
        errs[idx] = accel_l2_error(ref, bh, n);
        printf("  theta = %.1f  ->  error L2 = %.6e\n", thetas[idx], errs[idx]);
        octree_free(t);
        free(bh);
    }

    int monotone = 1;
    for (int idx = 1; idx < 4; idx++)
        if (errs[idx] > errs[idx-1]) monotone = 0;

    int exact = errs[3] < 1e-9;  /* theta=0 reproduce el metodo directo */

    printf("  Monotono decreciente: %s\n", monotone ? "si" : "no");
    printf("  theta=0 reproduce directo (err < 1e-9): %s\n", exact ? "si" : "no");

    int pass = monotone && exact;
    printf("  Resultado: %s\n\n", pass ? "PASS" : "FAIL");

    free(ref);
    return pass;
}

int run_all_tests(void) {
    printf("\n========== SUITE DE VALIDACION ==========\n\n");

    int total = 8, passed = 0;
    passed += test_two_body_orbit();
    passed += test_energy_conservation();
    passed += test_momentum_conservation();
    passed += test_morton_ordering();
    passed += test_tree_mass();
    passed += test_tree_cm();
    passed += test_bh_force_error();
    passed += test_bh_theta_convergence();

    printf("========== RESULTADO: %d/%d tests pasaron ==========\n\n", passed, total);
    return (passed == total) ? 0 : 1;
}
