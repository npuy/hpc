#include "mpi_types.h"

#ifdef USE_MPI

#include "particle.h"

void mpi_particle_type_init(MPI_Datatype *type) {
    /* Inicializado solo para no leer memoria indeterminada: MPI_Get_address
       usa la direccion de los campos, nunca su contenido. */
    Particle probe = {{0}, {0}, {0}, 0.0, 0, 0};

    /*
     * Tres bloques: los 10 doubles contiguos (pos[3], vel[3], acc[3], mass),
     * la clave Morton y el id. Los desplazamientos se miden sobre una instancia
     * real en vez de asumirlos, para que el padding que elija el compilador
     * quede descrito correctamente.
     */
    int          blocklen[3] = {10, 1, 1};
    MPI_Datatype types[3]    = {MPI_DOUBLE, MPI_UINT64_T, MPI_INT64_T};
    MPI_Aint     disp[3], base;

    MPI_Get_address(&probe,          &base);
    MPI_Get_address(&probe.pos[0],   &disp[0]);
    MPI_Get_address(&probe.morton,   &disp[1]);
    MPI_Get_address(&probe.id,       &disp[2]);
    for (int i = 0; i < 3; i++) disp[i] -= base;

    MPI_Datatype tmp;
    MPI_Type_create_struct(3, blocklen, disp, types, &tmp);
    /* Extent = sizeof(Particle): imprescindible para enviar arreglos. */
    MPI_Type_create_resized(tmp, 0, sizeof(Particle), type);
    MPI_Type_free(&tmp);
    MPI_Type_commit(type);
}

void mpi_particle_type_free(MPI_Datatype *type) {
    MPI_Type_free(type);
}

#endif /* USE_MPI */
