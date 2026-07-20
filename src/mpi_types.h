#ifndef MPI_TYPES_H
#define MPI_TYPES_H

/*
 * mpi_types.h — Tipo derivado de MPI para enviar structs Particle.
 *
 * Particle mezcla doubles con enteros de 64 bits, y el compilador puede
 * insertar padding entre los campos o al final del struct. Enviarlo como un
 * bloque de MPI_BYTE "funcionaría" solo si todos los nodos son binariamente
 * idénticos, y oculta errores de alineación en cuanto deja de serlo.
 *
 * En su lugar se construye un tipo derivado que describe explícitamente los
 * desplazamientos reales de cada campo (MPI_Type_create_struct) y se fija su
 * extent en sizeof(Particle) (MPI_Type_create_resized). Sin ese resize, MPI
 * calcula el extent a partir del último campo y el envío de un arreglo de
 * partículas quedaría desalineado a partir del segundo elemento.
 */

#ifdef USE_MPI

#include <mpi.h>

/*
 * Construye y hace commit del tipo derivado para Particle.
 * El caller debe liberarlo con mpi_particle_type_free antes de MPI_Finalize.
 */
void mpi_particle_type_init(MPI_Datatype *type);

/* Libera el tipo derivado creado por mpi_particle_type_init. */
void mpi_particle_type_free(MPI_Datatype *type);

#endif /* USE_MPI */

#endif
