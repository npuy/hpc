# Guía de ejecución en las máquinas pcunix

Todo lo de acá está **verificado en vivo**, no supuesto. Corrige varias cosas que
`mpi-setup.md` daba por buenas.

## Topología de acceso

```
portátil  --ssh-->  lulu (login-ens.fing.edu.uy)  --ssh-->  pcunixNN
```

`lulu` está en `~/.ssh/config` con su clave. **lulu es solo punto de acceso**: al
entrar muestra un banner que dice explícitamente que no está previsto ejecutar
programas ahí. Se usa como salto y nada más.

## Estado real del cluster

| Nodo | Estado |
|---|---|
| pcunix40 | **CAÍDA** (`No route to host`) — es la que `mpi-setup.md` usa de ejemplo |
| pcunix41 | OK |
| pcunix42 | responde, pero **conflicto de host key** (la clave ED25519 difiere de la de la IP) |
| pcunix43 | OK |
| pcunix44 | OK |
| pcunix45 | OK |

Por eso el script descubre los nodos vivos en tiempo de ejecución en vez de
fijarlos: cambian.

### Hardware por nodo

```
Intel(R) Core(TM) i3-4170 CPU @ 3.70GHz
2 núcleos FÍSICOS, 4 hilos (hyperthreading), 1 socket, 1 nodo NUMA
15 GB RAM
```

**Esto contradice el supuesto de tres planes seguidos.** Yo venía escribiendo que
había que "re-medir en pcunix con núcleos homogéneos" pensando en máquinas con 8 o
16 núcleos. Son 2 núcleos físicos por nodo: el portátil (4 performance + 4
efficiency) tiene más cómputo que un nodo pcunix.

**Lo que pcunix sí aporta, y el portátil no puede dar:**

1. **8 núcleos físicos homogéneos** (4 nodos × 2), sin la asimetría
   performance/efficiency que invalidó las tablas de las semanas 3-5.
2. **Ejecución realmente distribuida** entre máquinas, con la red de por medio.
   Es la deuda pendiente desde la semana 3.
3. **Portabilidad**: otro compilador, otro MPI, otra arquitectura.

Lo que **no** se va a poder medir acá es escalabilidad OpenMP más allá de 2 hilos
reales. Conviene decirlo en el informe en vez de presentar 4 hilos como si fueran
paralelismo.

### Software

```
gcc (GCC) 15.2.1 20260123 (Red Hat 15.2.1-7)
MPICH, vía módulo, en /usr/lib64/mpich/bin/mpicc
```

Distinto del portátil en las tres dimensiones que importan: compilador (clang →
gcc), implementación MPI (Open MPI → MPICH) y arquitectura (ARM → x86-64).

### El home es NFS compartido

Verificado escribiendo un archivo en lulu y leyéndolo desde pcunix41 y pcunix43.
Tiene tres consecuencias prácticas:

1. **Se compila una sola vez** y todos los nodos ven el binario. No hay que
   distribuir nada.
2. `~/.ssh/known_hosts` **también es compartido**, y por eso `mpirun` multi-nodo no
   pide confirmación de claves: se poblaron solas al hacer el barrido inicial
   desde lulu. Es la solución al problema que `mpi-setup.md` describe como
   "conectarse a mano a cada host y aceptar con yes".
3. El código se sube a `~/hpc/` desde el portátil y queda disponible en todo el
   cluster.

## Configuración obligatoria

```bash
module load mpi/mpich-x86_64
export FI_PROVIDER=tcp        # sin esto: "OFI endpoint open failed"
export OMP_PROC_BIND=close    # evita que el SO migre hilos y ensucie las medidas
export OMP_PLACES=cores
```

Nota: el `.bashrc` ya tiene el `module load`, y **falla en lulu** (`module: command
not found`) porque lulu no tiene módulos. Es ruido inofensivo, pero ensucia la
salida de todo comando por ssh.

## Colocación de procesos

MPICH reparte **round-robin** entre los hosts del hostfile. Verificado: con 4
hosts y `-np 8`, quedan exactamente 2 procesos por nodo, tanto con hostfile plano
como con sintaxis `host:2`.

```bash
printf 'pcunix41\npcunix43\npcunix44\npcunix45\n' > ~/hpc/hostfile
mpirun -np 8 -f ~/hpc/hostfile ./nbody_mpi ...
```

**Consecuencia para el diseño de los experimentos:** con 2 núcleos físicos por
nodo, correr `P > 2` en una sola máquina mide contención de hyperthreading — el
mismo defecto que invalidó las tablas del portátil. Los barridos de escalabilidad
tienen que ir **siempre distribuidos**, de modo que `P = 1..8` sean 1 a 8 núcleos
físicos, uno por proceso.

Es la razón por la que el script no usa nunca `mpirun -np 4` a secas.

## Puesta a punto

```bash
# desde el portátil, en la raíz del repo
tar czf /tmp/hpc-src.tgz src Makefile
scp /tmp/hpc-src.tgz lulu:~/
scp scripts/pcunix-bench.sh lulu:~/hpc/

ssh lulu
  rm -rf ~/hpc && mkdir -p ~/hpc
  tar xzf ~/hpc-src.tgz -C ~/hpc
  chmod +x ~/hpc/pcunix-bench.sh
  ssh pcunix41
    module load mpi/mpich-x86_64 && export FI_PROVIDER=tcp
    cd ~/hpc && make
```

El `tar` de macOS mete atributos extendidos que GNU tar reporta como
`Ignoring unknown extended header keyword 'LIBARCHIVE.xattr...'`. Es inofensivo.

## Qué ya se verificó

Compilación limpia, sin warnings, con GCC 15.2 + MPICH en x86-64:

| Configuración | Resultado |
|---|---|
| `./nbody --validate` | **9/9** |
| `mpirun -np 2 ./nbody_mpi --validate` | **19/19** |
| `mpirun -np 3 ./nbody_mpi --validate` | **19/19** |
| `mpirun -np 4 ./nbody_mpi --validate` | **19/19** |

Es un resultado de portabilidad en sí mismo: otro compilador, otro MPI y otra
arquitectura dan exactamente los mismos veredictos, incluidos los dos tests de
igualdad **bit a bit** (el 9, determinismo de OpenMP, y el 14, LET con θ=0).

P=3 importa aparte: verifica el reparto proporcional de ORB cuando el número de
procesos no es potencia de 2.

## Cómo sería la ejecución completa

`scripts/pcunix-bench.sh` está escrito y ya copiado a `~/hpc/`. **Todavía no se
ejecutó.** Descubre los nodos vivos, arma el hostfile, y corre ocho bloques:

| # | Bloque | Qué responde |
|---|---|---|
| 1 | Validación (secuencial, P=2/3/4/8 distribuido) | ¿sigue todo correcto en 8 núcleos y con la red de por medio? |
| 2 | ORB vs Morton, `n_ghost` vs N | **el resultado central de la semana 5**: exponente 0,66 vs 0,92 |
| 3 | Escalabilidad fuerte distribuida, P=1..8 | speedup real con núcleos homogéneos |
| 4 | Escalabilidad débil distribuida, N/P=12,5K | la que venía dando 5-7% en el portátil |
| 5 | OpenMP 1/2/4 hilos | con la salvedad de que 4 es hyperthreading |
| 6 | Híbrido 8×1 vs 4×2 | ¿conviene proceso o hilo? |
| 7 | Balance conteo vs costo | la inversión de los desbalances |
| 8 | Plummer vs uniform | verificación del diagnóstico de la semana 4 |

Cada corrida imprime una línea con `total`, `arbol`, `fuerzas`, `let`, `migr`,
`ghosts`, los dos desbalances y el `drift`. La salida completa queda en
`~/hpc/resultados-pcunix.txt`.

### Duración: no la sé, y hay que calibrarla

**No tengo una medición de cuánto tarda una corrida en un i3-4170**, así que
cualquier estimación que diera sería inventada. El único dato real: la validación
secuencial más P=2, 3 y 4 tardó más de dos minutos, y P=8 en un solo nodo (8
procesos sobre 2 núcleos físicos) no había terminado cuando lo corté.

Por eso conviene **correr primero un subconjunto corto** —bloques 1 y 2, que son
validación y el resultado central— cronometrarlo, y recién ahí decidir sobre la
matriz completa. Son máquinas compartidas: hay 3 usuarios conectados en cada una,
aunque con carga 0.00.

Si hace falta acortar, los tamaños del script (`N=50000` en los barridos,
`N=100000` en el bloque 2) son el parámetro a bajar primero.

### Ejecución

```bash
ssh lulu
  ssh pcunix41
    cd ~/hpc && ./pcunix-bench.sh              # en primer plano, se ve todo
    # o, para dejarlo corriendo:
    nohup ./pcunix-bench.sh > /dev/null 2>&1 &

# recuperar resultados desde el portátil
scp lulu:~/hpc/resultados-pcunix.txt docs/
```

## Qué mirar en los resultados

1. **19/19 en todas las configuraciones**, sobre todo P=8 distribuido. Si algo
   falla ahí y no en el portátil, es un problema real de la ejecución
   multi-máquina, no ruido.
2. **Bloque 2**: el ratio `ghosts/n_local` de ORB tiene que caer al crecer N,
   mientras el de Morton se queda clavado en ~2,5. Es el resultado de la semana 5
   y debería reproducirse en otra arquitectura.
3. **Bloque 3**: `arbol` con ORB debe mantenerse plano al crecer P; con Morton
   debe crecer. Es el criterio que la semana 4 no cumplía.
4. **Bloque 4**: la eficiencia débil venía dando 5-7% en el portátil. Con núcleos
   homogéneos y sin sobresuscribir debería mejorar; **si no mejora, el problema es
   algorítmico y no del hardware**, y eso es un hallazgo que hay que reportar.
5. **Bloque 5**: esperar ~1,9× con 2 hilos y poca ganancia con 4. Si 4 hilos
   dieran mucho más que 2, algo está mal medido.

## Correcciones a `mpi-setup.md`

| Dice | Realidad medida |
|---|---|
| Ejemplos con `pcunix40,pcunix42` | pcunix40 está caída; pcunix42 tiene conflicto de host key. Usar 41, 43, 44, 45 |
| Hay que aceptar claves SSH a mano en cada host | El `known_hosts` es NFS compartido: se pueblan solas desde lulu |
| (no lo menciona) | El home es compartido: se compila una sola vez |
| (no lo menciona) | 2 núcleos físicos por nodo; `P>2` en un nodo sobresuscribe |
| (no lo menciona) | lulu es solo punto de acceso, no se ejecuta ahí |
