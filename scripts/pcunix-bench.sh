#!/bin/bash
# Benchmarks en las maquinas pcunix de fing.
#
# Entorno (ver docs/mpi-setup.md): MPICH por modulo, FI_PROVIDER=tcp obligatorio,
# home NFS compartido entre lulu y todas las pcunix -- se compila una sola vez y
# todos los nodos ven el binario, y ~/.ssh/known_hosts tambien es compartido, que
# es por lo que mpirun multi-nodo no pide confirmacion de claves.
#
# HARDWARE, que condiciona como se leen los resultados: cada pcunix es un Intel
# i3-4170 con 2 nucleos FISICOS y 4 hilos por hyperthreading. Poco computo por
# nodo, pero 4 nodos dan 8 nucleos fisicos HOMOGENEOS, que es justo lo que el
# portatil (4 performance + 4 efficiency) no puede dar.
#
# De ahi la decision de colocacion: los barridos de escalabilidad van SIEMPRE
# distribuidos, 2 procesos por nodo (MPICH reparte round-robin), de modo que
# P=1..8 son 1..8 nucleos fisicos sin sobresuscribir. Correr P>2 en un solo nodo
# mide contencion de hyperthreading, que es el mismo defecto que invalidaba las
# tablas del portatil.
#
# Uso:
#   ./pcunix-bench.sh                  todos los bloques
#   ./pcunix-bench.sh 1,2              solo los bloques 1 y 2 (calibracion)
#   ./pcunix-bench.sh 3-5              rango de bloques
#   ./pcunix-bench.sh 1,2 salida.txt   eligiendo el archivo de salida
#
# Bloques: 1 validacion | 2 ORB vs Morton | 3 escalab. fuerte | 4 escalab. debil
#          5 OpenMP | 6 hibrido | 7 balance | 8 plummer vs uniform

set -u

# Expande "1,2" o "3-5" a una lista de bloques a correr. Sin argumento, todos.
BLOQUES="${1:-1-8}"
case "$BLOQUES" in
    *[!0-9,-]*) BLOQUES="1-8"; OUT="${1:-}" ;;   # el 1er arg era el archivo
esac
OUT="${2:-${OUT:-$HOME/hpc/resultados-pcunix.txt}}"
[ -z "$OUT" ] && OUT="$HOME/hpc/resultados-pcunix.txt"

SELECTED=""
IFS=, read -ra _parts <<< "$BLOQUES"
for _p in "${_parts[@]}"; do
    case "$_p" in
        *-*) for _i in $(seq "${_p%%-*}" "${_p##*-}"); do SELECTED="$SELECTED $_i"; done ;;
        *)   SELECTED="$SELECTED $_p" ;;
    esac
done
# quiere <n> -> 0 si el bloque n esta seleccionado
quiere() { case " $SELECTED " in *" $1 "*) return 0;; *) return 1;; esac; }

module load mpi/mpich-x86_64 2>/dev/null
export FI_PROVIDER=tcp
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=1

cd "$HOME/hpc" || exit 1

# Nodos vivos: se descubren en vivo y no se fijan. pcunix40, la que usa como
# ejemplo docs/mpi-setup.md, esta caida.
NODES=""
for h in pcunix41 pcunix42 pcunix43 pcunix44 pcunix45; do
    if timeout 6 ssh -o BatchMode=yes -o StrictHostKeyChecking=no "$h" true 2>/dev/null; then
        NODES="$NODES $h"
    fi
done
NODES=$(echo $NODES)
NNODES=$(echo $NODES | wc -w)
HF="$HOME/hpc/hostfile"
printf '%s\n' $NODES > "$HF"
MAXP=$((NNODES * 2))     # 2 procesos por nodo = 1 por nucleo fisico

exec > >(tee "$OUT") 2>&1

echo "=================================================================="
echo " Benchmarks pcunix  --  $(date)"
echo "=================================================================="
echo "Bloques a correr:$SELECTED"
echo "Nodos vivos ($NNODES): $NODES"
echo "Por nodo: $(nproc) hilos logicos, $(lscpu | grep -m1 '^Core(s)' | tr -dc 0-9) nucleos fisicos"
echo "Tope sin sobresuscribir: P=$MAXP"
lscpu | grep -m1 'Model name' | sed 's/  */ /g'
gcc --version | head -1
echo "MPI: $(mpirun --version 2>&1 | head -1)"
echo

# Extrae metricas de una corrida.
metrics() {
    awk '
        /^  arbol/                {tr=$4}
        /^  fuerzas/              {fo=$4}
        /^  let\(/                {le=$4}
        /^  migracion/            {mi=$4}
        /Tiempo de ejecucion/     {to=$4}
        /Fantasmas por proceso/   {split($6,a,"="); g=a[2]}
        /Desbalance de particulas/{dp=$5}
        /Desbalance de trabajo/   {dw=$5}
        /Energy drift/            {dr=$3}
        END{printf "total=%-8s arbol=%-8s fuerzas=%-8s let=%-8s migr=%-8s ghosts=%-8s desbP=%-7s desbW=%-7s drift=%s\n",
                   to,tr,fo,(le==""?"-":le),mi,(g==""?"-":g),dp,dw,dr}'
}

# run <nprocs> <args...>   -- siempre distribuido, 2 por nodo
run() { local p=$1; shift; mpirun -np "$p" -f "$HF" ./nbody_mpi "$@" --metrics 2>&1 | metrics; }

if quiere 1; then
echo "########## 1. Validacion ##########"
printf "secuencial      : "; ./nbody --validate 2>&1 | grep -oE "[0-9]+/[0-9]+ tests pasaron"
for P in 2 4 8; do
    [ "$P" -gt "$MAXP" ] && continue
    printf "MPI P=%-2s distrib: " $P
    mpirun -np $P -f "$HF" ./nbody_mpi --validate 2>&1 | grep -oE "[0-9]+/[0-9]+ tests pasaron"
done
printf "MPI P=3 distrib : "; mpirun -np 3 -f "$HF" ./nbody_mpi --validate 2>&1 | grep -oE "[0-9]+/[0-9]+ tests pasaron"
echo

fi

if quiere 2; then
echo "########## 2. ORB vs Morton: volumen del LET (P=4 distribuido) ##########"
echo "El resultado central de la semana 5: n_ghost ~ N^k, k=0.66 (orb) vs 0.92 (morton)."
for N in 12500 25000 50000 100000; do
    for D in morton orb; do
        printf "N=%-7s %-7s : " $N $D
        run 4 -n $N -t 0.003 --decomp $D
    done
done
echo

fi

if quiere 3; then
echo "########## 3. Escalabilidad fuerte distribuida (N=50K, 1 hilo) ##########"
echo "P=1..8 son 1..8 nucleos fisicos homogeneos repartidos en $NNODES nodos."
for P in 1 2 4 8; do
    [ "$P" -gt "$MAXP" ] && continue
    for D in morton orb; do
        printf "P=%-2s %-7s : " $P $D
        run $P -n 50000 -t 0.005 --decomp $D
    done
done
echo

fi

if quiere 4; then
echo "########## 4. Escalabilidad debil distribuida (N/P = 12.5K) ##########"
for P in 1 2 4 8; do
    [ "$P" -gt "$MAXP" ] && continue
    for D in morton orb; do
        printf "P=%-2s N=%-7s %-7s : " $P $((12500*P)) $D
        run $P -n $((12500*P)) -t 0.005 --decomp $D
    done
done
echo

fi

if quiere 5; then
echo "########## 5. Escalabilidad OpenMP (1 proceso, 1 nodo, N=50K) ##########"
echo "OJO: 2 nucleos fisicos. 4 hilos es hyperthreading, no paralelismo real."
for T in 1 2 4; do
    printf "hilos=%-2s : " $T
    OMP_NUM_THREADS=$T mpirun -np 1 ./nbody_mpi -n 50000 -t 0.005 --metrics 2>&1 | metrics
done
export OMP_NUM_THREADS=1
echo

fi

if quiere 6; then
echo "########## 6. Hibrido a 8 nucleos fisicos (N=50K) ##########"
for cfg in "8 1" "4 2"; do
    set -- $cfg
    [ "$1" -gt "$MAXP" ] && continue
    printf "%s proc x %s hilos : " $1 $2
    OMP_NUM_THREADS=$2 mpirun -np $1 -f "$HF" ./nbody_mpi -n 50000 -t 0.005 --metrics 2>&1 | metrics
done
export OMP_NUM_THREADS=1
echo

fi

if quiere 7; then
echo "########## 7. Balance: conteo vs costo (P=8 distribuido, N=50K) ##########"
for B in count work; do
    printf "%-6s : " $B
    run $((MAXP < 8 ? MAXP : 8)) -n 50000 -t 0.01 --balance $B
done
echo

fi

if quiere 8; then
echo "########## 8. Plummer vs uniform (P=4 distribuido, N=50K) ##########"
for I in plummer uniform; do
    for D in morton orb; do
        printf "%-8s %-7s : " $I $D
        run 4 -n 50000 -t 0.003 --init $I --decomp $D
    done
done
echo

fi

echo "=================================================================="
echo " Fin  --  $(date)"
echo "=================================================================="
