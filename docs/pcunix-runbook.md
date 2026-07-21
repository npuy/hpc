# Runbook: ejecutar las pruebas en pcunix, paso a paso

Guía operativa. Cada paso trae el comando exacto, qué se espera ver, y qué hacer
si falla. Pensada para seguirse desde cero sin contexto previo.

Para entender **por qué** el cluster se usa así, ver
[pcunix-guia.md](pcunix-guia.md). Para entender qué hace el programa, ver
[arquitectura.md](arquitectura.md).

---

## Paso 0 — Antes de empezar

Hace falta:

- Acceso a `lulu` configurado en `~/.ssh/config` (ya está).
- El repo local en `~/Documents/facu/hpc`.

Dos cosas que conviene saber de entrada, porque si no parecen errores:

- **Todo comando por ssh imprime un banner** de varias líneas y un
  `module: command not found`. Es ruido normal: `lulu` no tiene módulos y el
  `.bashrc` intenta cargarlos igual.
- **`lulu` es solo punto de acceso.** No se compila ni se ejecuta nada ahí.

---

## Paso 1 — Verificar acceso

```bash
ssh lulu 'hostname'
```

**Esperado:** entre el banner, aparece `login-ens.fing.edu.uy`.

**Si falla:** revisar que exista `~/.ssh/fing-lulu/id_ed25519` y que el host
`lulu` esté en `~/.ssh/config`.

---

## Paso 2 — Ver qué nodos están vivos

Las máquinas se caen y vuelven; no asumir cuáles hay.

```bash
ssh lulu 'for h in pcunix41 pcunix42 pcunix43 pcunix44 pcunix45; do
  printf "%-10s " $h
  ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=5 $h \
      "echo OK \$(nproc)hilos carga=\$(cut -d\" \" -f1 /proc/loadavg)" 2>&1 | tail -1
done'
```

**Esperado:** al menos 3 nodos con `OK 4hilos carga=0.00`.

**Notas:**
- `pcunix40` está caída (`No route to host`). No está en la lista a propósito.
- `pcunix42` suele dar conflicto de host key; si aparece, ignorarla y seguir con
  las demás. Con 3 nodos alcanza.
- Este paso además **acepta las claves SSH** de cada nodo. Como el `known_hosts`
  vive en el home NFS compartido, queda aceptado para todos los nodos a la vez, y
  es lo que hace que `mpirun` multi-máquina no pida confirmación después.
- Si la carga no es 0.00, hay alguien más usando la máquina: las mediciones de
  tiempo van a salir sucias. Conviene esperar o elegir otros nodos.

---

## Paso 3 — Subir el código

Desde la raíz del repo, en el portátil:

```bash
cd ~/Documents/facu/hpc
tar czf /tmp/hpc-src.tgz src Makefile
scp /tmp/hpc-src.tgz lulu:~/
scp scripts/pcunix-bench.sh lulu:~/

ssh lulu 'rm -rf ~/hpc && mkdir -p ~/hpc &&
          tar xzf ~/hpc-src.tgz -C ~/hpc &&
          mv ~/pcunix-bench.sh ~/hpc/ && chmod +x ~/hpc/pcunix-bench.sh &&
          rm ~/hpc-src.tgz && ls ~/hpc'
```

**Esperado:** `Makefile  pcunix-bench.sh  src`

**Nota:** el `tar` de macOS mete atributos extendidos y GNU tar avisa
`Ignoring unknown extended header keyword 'LIBARCHIVE.xattr...'` una vez por
archivo. Es inofensivo.

Como el home es NFS compartido, con subirlo una vez queda visible en todos los
nodos. No hay que copiar nada a cada máquina.

---

## Paso 4 — Compilar

```bash
ssh lulu 'ssh pcunix41 "bash -lc \"
  module load mpi/mpich-x86_64
  cd ~/hpc && make clean >/dev/null 2>&1
  make 2>&1 | grep -Ei \\\"error|warning\\\"
  ls -la nbody nbody_mpi
\""'
```

**Esperado:** ningún error ni warning, y los dos binarios listados.

**Si falla `mpicc: command not found`:** faltó el `module load`. El `bash -lc`
está para que se lea el `.bashrc`, que ya lo tiene.

---

## Paso 5 — Prueba de humo (2 minutos)

Antes de nada largo, confirmar que corre:

```bash
ssh lulu 'ssh pcunix41 "bash -lc \"
  module load mpi/mpich-x86_64; export FI_PROVIDER=tcp; cd ~/hpc
  mpirun -np 2 ./nbody_mpi -n 10000 -t 0.005 --metrics 2>&1 |
    grep -E \\\"Particulas:|Checksum|Error de masa|drift\\\"
\""'
```

**Esperado:**

```
Energy drift: ~1e-05
Particulas: 10000 / 10000  (conservadas)
Checksum id: OK
Error de masa del LET: ~1e-14
```

**Si aparece `OFI endpoint open failed ... Cannot allocate memory`:** faltó
`export FI_PROVIDER=tcp`. Es obligatorio y no es opcional en estas máquinas.

---

## Paso 6 — Calibrar cuánto tarda

**No sé cuánto tarda una corrida en un i3-4170 y no quiero inventarlo.** Este paso
lo mide, para decidir con un número real si conviene lanzar la matriz completa.

```bash
ssh lulu 'ssh pcunix41 "bash -lc \"
  module load mpi/mpich-x86_64; export FI_PROVIDER=tcp; cd ~/hpc
  printf \\\"pcunix41\\npcunix43\\npcunix44\\npcunix45\\n\\\" > hostfile
  echo UNA CORRIDA TIPICA:
  time mpirun -np 4 -f hostfile ./nbody_mpi -n 50000 -t 0.005 > /dev/null
\""'
```

**Cómo usarlo:** el script completo hace **unas 60 corridas** de tamaño parecido o
menor. Multiplicar el `real` que salga por 60 da el orden de magnitud del total.

| Si una corrida tarda | El script completo tarda | Recomendación |
|---|---|---|
| < 5 s | < 5 min | lanzar todo (Paso 8) |
| 5–20 s | 5–20 min | lanzar todo, en background |
| > 20 s | > 20 min | bajar `N` en el script, o correr por bloques |

Para bajar los tamaños: editar `~/hpc/pcunix-bench.sh` y cambiar los `-n 50000`
por `-n 25000` y el `100000` del bloque 2 por `50000`.

---

## Paso 7 — Corrida corta: validación y el resultado central

Antes de la matriz completa, los dos bloques que más importan.

```bash
ssh lulu 'ssh pcunix41 "bash -lc \"cd ~/hpc && ./pcunix-bench.sh 1,2\""'
```

El script acepta bloques sueltos: `1,2` o rangos `3-5` o mezclas `1,3-4,8`.

**Bloque 1 — Validación.** Esperado: `9/9` secuencial y `19/19` en todas las
configuraciones MPI, incluida P=8 distribuida en 4 máquinas.

> Si algo falla acá y **no** falla en el portátil, es un problema real de la
> ejecución multi-máquina, no ruido. Vale la pena parar y mirarlo.

**Bloque 2 — ORB vs Morton.** Esperado: la columna `ghosts` de `orb` crece mucho
más despacio que la de `morton` al subir N. En el portátil, a P=4:

| N | morton | orb |
|---|---|---|
| 12.500 | 8.324 | 4.281 |
| 25.000 | 16.271 | 6.830 |
| 50.000 | 33.453 | 10.789 |
| 100.000 | 61.056 | 16.931 |

Si se reproduce la forma (morton ~ lineal, orb claramente sublineal), el
resultado de la semana 5 vale también en otra arquitectura.

---

## Paso 8 — Corrida completa

Con la calibración hecha, lanzar todo. En background, para que no dependa de que
la sesión ssh siga viva:

```bash
ssh lulu 'ssh pcunix41 "bash -lc \"
  cd ~/hpc && nohup ./pcunix-bench.sh > /dev/null 2>&1 & echo lanzado
\""'
```

**Ver el avance:**

```bash
ssh lulu 'tail -20 ~/hpc/resultados-pcunix.txt'
```

**Ver si todavía corre:**

```bash
ssh lulu 'ssh pcunix41 "pgrep -u \$(whoami) -f nbody_mpi >/dev/null && echo CORRIENDO || echo TERMINADO"'
```

**Terminó cuando** la última línea del archivo dice `Fin -- <fecha>`.

---

## Paso 9 — Traer los resultados

```bash
cd ~/Documents/facu/hpc
scp lulu:~/hpc/resultados-pcunix.txt docs/
```

Quedan en `docs/resultados-pcunix.txt`, listos para armar las tablas del informe.

---

## Paso 10 — Cortar una corrida

Si algo se cuelga o hay que liberar las máquinas:

```bash
ssh lulu 'for h in pcunix41 pcunix43 pcunix44 pcunix45; do
  ssh $h "pkill -u \$(whoami) -f nbody_mpi; pkill -u \$(whoami) -f pcunix-bench" 2>/dev/null
done; echo limpiado'
```

Matar los procesos en **todos** los nodos, no solo en el que lanzó: `mpirun` deja
procesos en cada máquina.

---

## Qué mirar en los resultados

| Bloque | Qué debería verse | Qué significaría lo contrario |
|---|---|---|
| 1 | 19/19 en todo | Falla propia de multi-máquina |
| 2 | `ghosts` de orb crece sublineal | El resultado de la semana 5 no es portable |
| 3 | `arbol` con orb plano al subir P; con morton, creciendo | — |
| 4 | Eficiencia débil mejor que el 5-7% del portátil | **Si no mejora, el problema es algorítmico y no del hardware. Es un hallazgo, no un fracaso.** |
| 5 | ~1,9× con 2 hilos, poca ganancia con 4 | Si 4 hilos dieran mucho más que 2, está mal medido: hay 2 núcleos físicos |
| 7 | Desbalance de trabajo baja, el de partículas sube | La inversión es la señal de que el balance por costo funciona |

---

## Problemas conocidos

| Síntoma | Causa | Solución |
|---|---|---|
| `OFI endpoint open failed ... Cannot allocate memory` | Falta forzar TCP | `export FI_PROVIDER=tcp` |
| `mpicc: command not found` | Falta el módulo | `module load mpi/mpich-x86_64`, y usar `bash -lc` |
| `No route to host` a pcunix40 | Está caída | Usar 41, 43, 44, 45 |
| `host key ... differs from the key for the IP` en pcunix42 | Clave inconsistente | Saltearla; con 3 nodos alcanza |
| `Host key verification failed` en `mpirun` | Claves no aceptadas | Correr el Paso 2, que las acepta todas de una |
| `Connection timed out during banner exchange` | Se intentó `ssh -J lulu pcunixNN` contra un nodo caído | Usar ssh anidado (`ssh lulu` y de ahí `ssh pcunixNN`). *No verifiqué si `-J` funciona contra un nodo vivo; el único intento fue contra pcunix40, que está caída.* |
| `module: command not found` al entrar a lulu | lulu no tiene módulos | Ruido inofensivo, ignorar |
| `Ignoring unknown extended header keyword` en el tar | Atributos extendidos de macOS | Ruido inofensivo, ignorar |
| Tiempos erráticos entre corridas | Otro usuario en la máquina | Verificar carga (Paso 2) y repetir |

---

## Referencia rápida de comandos

```bash
# subir cambios de código y recompilar
cd ~/Documents/facu/hpc && tar czf /tmp/hpc-src.tgz src Makefile
scp /tmp/hpc-src.tgz lulu:~/ && ssh lulu 'tar xzf ~/hpc-src.tgz -C ~/hpc && rm ~/hpc-src.tgz'
ssh lulu 'ssh pcunix41 "bash -lc \"cd ~/hpc && make\""'

# correr bloques sueltos
ssh lulu 'ssh pcunix41 "bash -lc \"cd ~/hpc && ./pcunix-bench.sh 1,2\""'

# una corrida puntual
ssh lulu 'ssh pcunix41 "bash -lc \"
  module load mpi/mpich-x86_64; export FI_PROVIDER=tcp; cd ~/hpc
  mpirun -np 8 -f hostfile ./nbody_mpi -n 50000 -t 0.005 --metrics\""'

# traer resultados
scp lulu:~/hpc/resultados-pcunix.txt docs/
```
