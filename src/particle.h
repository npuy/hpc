#ifndef PARTICLE_H
#define PARTICLE_H

/*
 * particle.h — Estructura de datos central y funciones de inicialización/E/S.
 *
 * Define el tipo Particle (96 bytes) y provee generadores de condiciones
 * iniciales para distintos escenarios: sistema de 2 cuerpos (validación),
 * esfera de Plummer (simulación astrofísica realista) y cubo uniforme
 * (benchmark simple). Incluye lectura/escritura en CSV para snapshots.
 *
 * El campo id es la identidad global estable de la partícula: sobrevive al
 * reordenamiento por Morton y a la migración entre procesos MPI, y es lo que
 * permite verificar que ninguna partícula se pierde ni se duplica (tests 11 y
 * 13). Los generadores lo asignan como el índice inicial 0..n-1.
 *
 * El campo work es el costo medido del recorrido Barnes-Hut de la partícula en
 * el último paso (cantidad de interacciones evaluadas). Vive dentro del struct
 * y no en un arreglo paralelo a propósito: así viaja solo en la migración, sin
 * que haya que mantener sincronizada una segunda permutación. Es el peso que
 * usa el rebalanceo dinámico (ver domain.h).
 *
 * Convención de la semana 4: una partícula con id == -1 es un FANTASMA, es
 * decir un nodo remoto importado por el LET (ver let.h). Los fantasmas ejercen
 * fuerza pero no la reciben, y quedan excluidos de checksums, energías y
 * migración.
 */

#include <stdint.h>

#define PARTICLE_GHOST_ID ((int64_t)-1)

typedef struct {
    double pos[3];    /* posición en espacio 3D (x, y, z) */
    double vel[3];    /* velocidad (vx, vy, vz) */
    double acc[3];    /* aceleración gravitatoria acumulada (ax, ay, az) */
    double mass;      /* masa de la partícula (masa total del sistema = 1) */
    uint64_t morton;  /* clave Morton de 63 bits para ordenamiento Z-order */
    int64_t  id;      /* identidad global estable; -1 marca un fantasma del LET */
    int64_t  work;    /* interacciones del último paso: peso del rebalanceo */
} Particle;

/*
 * Configura un sistema de 2 cuerpos de masa M/2 en órbita circular simétrica
 * respecto al origen. Separación = 1, velocidades calculadas analíticamente
 * para G=1. Escribe n=2 en *n. Útil para validar el integrador.
 */
void init_two_body(Particle *p, int *n);

/*
 * Genera n partículas según el modelo de Plummer: posiciones muestreadas por
 * CDF inversa, velocidades por rejection sampling (Aarseth et al. 1974).
 * Masa total = 1, radio de escala a = 1. Aplica corrección de centro de masa
 * para centrar el sistema en el origen con momento total nulo.
 */
void init_plummer(Particle *p, int n, unsigned seed);

/*
 * Genera n partículas con posiciones uniformes en el cubo [-0.5, 0.5]³ y
 * velocidades nulas. Masa total = 1 repartida equitativamente. Útil como
 * benchmark de colapso gravitacional (cold collapse).
 */
void init_uniform_cube(Particle *p, int n, unsigned seed);

/*
 * Serializa el estado de n partículas al archivo fname en CSV con header
 * "t,id,mass,x,y,z,vx,vy,vz" y precisión de 12 dígitos significativos.
 * El parámetro t marca el instante temporal del snapshot.
 */
void write_particles_csv(const char *fname, const Particle *p, int n, double t);

/*
 * Deserializa partículas desde un archivo CSV con el formato de
 * write_particles_csv. Salta la línea de header. Escribe la cantidad
 * leída en *n. Retorna 0 en éxito, -1 si no puede abrir el archivo.
 * El caller debe asegurar que p tenga espacio suficiente.
 */
int  read_particles_csv(const char *fname, Particle *p, int *n);

#endif
