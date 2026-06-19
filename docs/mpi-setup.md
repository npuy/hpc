# Configuracion y ejecucion de MPI en maquinas de fing

## Entorno de ejecucion

Las maquinas de fing (pcunix40, pcunix42, etc.) usan modulos para gestionar herramientas.

## Configuracion de MPI

```bash
# Cargar modulo (una vez por sesion)
module load mpi/mpich-x86_64

# Para carga automatica en cada login
echo "module load mpi/mpich-x86_64" >> ~/.bashrc
```

**Fix obligatorio**: MPICH falla con el proveedor de red por defecto (PSM3/OFI). Hay que forzar TCP:

```bash
export FI_PROVIDER=tcp
# Agregar al .bashrc para que sea permanente
```

Sin esto se obtiene el error: `OFI endpoint open failed ... Cannot allocate memory`

## Compilacion

```bash
mpicc programa.c -o programa
```

## Ejecucion

```bash
# Local (1 maquina, 4 procesos)
mpirun -np 4 ./programa

# Distribuido (multiples maquinas)
mpirun -np 8 -hosts pcunix40,pcunix42 ./programa

# Distribuido con archivo de hosts (un host por linea)
mpirun -np 16 -hostfile mis_hosts ./programa
```

## Ejecucion distribuida - prerequisitos SSH

Antes de ejecutar en multiples maquinas, hay que aceptar las claves SSH de cada host. Desde cada maquina, conectarse a todas las demas (usando nombre corto y FQDN):

```bash
ssh pcunix42
ssh pcunix42.fing.edu.uy
```

Aceptar con `yes` en cada caso. Tambien verificar la conexion inversa (desde pcunix42 hacia pcunix40).

Si no se hace esto, `mpirun` falla con `Host key verification failed`.
