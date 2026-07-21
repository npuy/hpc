#include "mpi_types.h"

#ifdef USE_MPI

#include "particle.h"

void mpi_particle_type_init(MPI_Datatype *type) {
    /* Inicializado solo para no leer memoria indeterminada: MPI_Get_address
       usa la direccion de los campos, nunca su contenido. */
    Particle probe = {{0}, {0}, {0}, 0.0, 0, 0, 0};

    /*
     * Cuatro bloques: los 10 doubles contiguos (pos[3], vel[3], acc[3], mass),
     * la clave Morton, el id y el work. Los desplazamientos se miden sobre una
     * instancia real en vez de asumirlos, para que el padding que elija el
     * compilador quede descrito correctamente.
     *
     * id y work son ambos int64_t y con toda probabilidad quedan contiguos, pero
     * se describen como bloques separados a propósito: un bloque de 2 asumiría
     * que no hay padding entre ellos, y el costo de la suposición sería que la
     * migración copiara basura en silencio. Agregar un campo al struct y olvidar
     * este archivo es el error más fácil de cometer en la semana 4.
     */
    int          blocklen[4] = {10, 1, 1, 1};
    MPI_Datatype types[4]    = {MPI_DOUBLE, MPI_UINT64_T, MPI_INT64_T, MPI_INT64_T};
    MPI_Aint     disp[4], base;

    MPI_Get_address(&probe,          &base);
    MPI_Get_address(&probe.pos[0],   &disp[0]);
    MPI_Get_address(&probe.morton,   &disp[1]);
    MPI_Get_address(&probe.id,       &disp[2]);
    MPI_Get_address(&probe.work,     &disp[3]);
    for (int i = 0; i < 4; i++) disp[i] -= base;

    MPI_Datatype tmp;
    MPI_Type_create_struct(4, blocklen, disp, types, &tmp);
    /* Extent = sizeof(Particle): imprescindible para enviar arreglos. */
    MPI_Type_create_resized(tmp, 0, sizeof(Particle), type);
    MPI_Type_free(&tmp);
    MPI_Type_commit(type);
}

void mpi_particle_type_free(MPI_Datatype *type) {
    MPI_Type_free(type);
}

#endif /* USE_MPI */
