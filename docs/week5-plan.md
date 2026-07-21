# Semana 5 — Plan Detallado: ORB para que el LET comprima

## Objetivo

Reemplazar la descomposición por rangos de clave Morton por **ORB (bisección
recursiva ortogonal)**, de modo que cada proceso sea dueño de una **caja** en vez
de un tramo de curva. Con dominios que son cajas, el criterio de exportación del
LET pasa a ser exacto y ajustado, y el LET debería recuperar el comportamiento que
la semana 4 no consiguió.

Es una semana de **una sola corrección**, dirigida al resultado negativo medido:

| Criterio de la semana 4 | Umbral | Medido |
|---|---|---|
| `n_ghost` crece sublinealmente con N | k < 0.9 | **k = 0.92** |
| `tree_time` deja de crecer con P | no crece | **4,6×** de P=1 a P=8 |

## El diagnóstico que hay que atacar

De [week4-report.md](week4-report.md), medido y no inferido:

1. **La caja del proceso 3 cubre el 44,6% del volumen del dominio** (N=100K, P=4).
   Contra un destino así, `d_min ≈ 0` para casi todo nodo, el MAC nunca se cumple
   y se exporta todo sin resumir.
2. **`--theta-let` no tiene ningún efecto** (de 0,5 a 4,0 el volumen no se mueve):
   la rama que acepta el resumen prácticamente nunca se toma.
3. **Con `--init uniform` el volumen cae de 61.056 a 22.029 fantasmas**: la
   patología es de la distribución centralmente concentrada, no del código.

La causa raíz: en Plummer el pico de densidad cae justo donde se tocan los
octantes de nivel 1, así que las claves Morton del core se reparten por todo el
espacio de claves. Cualquier corte por percentil de conteo le da a cada proceso
una mezcla de core y halo, y **los dominios se interpenetran geométricamente**.

> **Por qué esto sí lo arregla ORB.** Con Morton, el dominio del proceso 3
> *contenía* partículas pegadas al core, así que su AABB contenía el core y
> ningún nodo del core podía resumirse para él. Con ORB los dominios son regiones
> disjuntas y contiguas: la caja de un proceso periférico **no contiene el core**,
> su cara interna está a una distancia real, y los nodos del core sí satisfacen el
> MAC contra ella. El tráfico que queda es de superficie, que es lo que un LET
> debe costar.

### Alternativas descartadas, y por qué

| Alternativa | Por qué no |
|---|---|
| Curva de Hilbert en vez de Morton | Mejor localidad que Morton, pero **un tramo contiguo sigue sin ser una región compacta**. Ataca el síntoma, no la causa. |
| Descriptor multi-caja (top-tree) del destino | Describiría mejor la región, pero si el dominio del proceso 3 genuinamente contiene partículas del core, el detalle del core hay que mandárselo igual. No cambia la conclusión. |
| Dual-tree traversal | Misma objeción que el anterior, con mucho más código. |

Ninguna de las tres toca la interpenetración de los dominios, que es el problema.

## Estado de partida

```
src/
  domain.h / .c       -- rangos Morton + biseccion por histograma   [REESCRIBIR]
  migration.h / .c    -- dueno por clave (busqueda binaria)         [MODIFICAR]
  let.h / .c          -- MAC conservador vs AABB del destino        [casi sin cambios]
  morton.h / .c       -- claves y orden                             [cambia de rol]
  octree.h / .c       -- octree_build_box, insert_particles         [sin cambios]
  validation.c        -- tests 12 y 16                              [MODIFICAR]
```

**Lo que NO se toca:** `let.c` (el criterio ya es correcto, lo que fallaba era la
geometría que recibía), `octree.c`, `barnes_hut.c`, el balance por costo.

---

## Día 1: El árbol ORB y la medición que decide

### Estudiante A — Particionamiento ORB

```c
typedef struct {
    double min[3], max[3];  /* caja de este nodo */
    int    first, nprocs;   /* tramo de ranks que cubre */
    int    axis;            /* eje de corte; -1 si es hoja */
    double split;           /* coordenada del corte */
    int    left, right;     /* indices de hijos; -1 si hoja */
} OrbNode;
```

El árbol tiene `2P-1` nodos y está **replicado en todos los procesos** — son unos
pocos cientos de bytes, así que no hay razón para distribuirlo.

Algoritmo, recursivo sobre el grupo de procesos:

```
orb_split(caja, [first, first+nprocs)):
    si nprocs == 1: hoja, este rank es dueno de la caja
    axis = dimension mas larga de la caja
    nL = nprocs / 2
    objetivo = trabajo_total * nL / nprocs
    bisecar c en [caja.min[axis], caja.max[axis]]:
        w_local  = Σ work de las particulas locales con pos[axis] < c
        MPI_Allreduce(w_local -> w_global)
        ajustar el intervalo segun w_global vs objetivo
    recursion en las dos mitades con las cajas partidas por c
```

Dos decisiones a documentar:

- **Eje = dimensión más larga.** Es la heurística estándar y evita cajas con forma
  de lámina, que tienen mucha superficie y por lo tanto mucho LET.
- **Se bisecta el trabajo, no el conteo.** El campo `work` de la semana 4 ya
  existe y ya está validado (test 17): ORB hereda el balance por costo sin código
  nuevo. Con `w_global == 0` (primer paso) se cae al conteo, igual que
  `domain_partition_ex`.

Detalles finos:

- **P no potencia de 2:** partir en `nL = nprocs/2` y `nprocs - nL`, con objetivo
  proporcional `w_total · nL / nprocs`. Sale gratis y evita una restricción fea.
- **El conteo por candidato no puede ser búsqueda binaria** — el arreglo local no
  está ordenado por la coordenada que se bisecta. Empezar con un **barrido lineal**
  (~40 iteraciones × n_local por nivel × log P niveles). Con n=25K y P=8 son ~3M
  operaciones por repartición, y las reparticiones son cada 20 pasos: no debería
  aparecer en las métricas. **Medirlo antes de optimizarlo**; si aparece, el
  reemplazo es un histograma de B bins con un solo Allreduce por nivel.
- **Terminación de la bisección:** todos los procesos derivan `lo`/`hi` solo de
  valores globales, así que recorren las mismas iteraciones y salen a la vez — el
  mismo argumento que ya hace válido a `domain_partition`. Cortar por número fijo
  de iteraciones (≈50), no por tolerancia relativa.
- **Partículas exactamente sobre el plano:** usar `<` de forma consistente en el
  conteo y en la asignación de dueño. Una discrepancia acá pierde partículas.

```c
void domain_partition_orb(Domain *d, const Particle *p, int n_local,
                          int64_t n_global, int use_work, MPI_Comm comm);
int  domain_owner_pos(const Domain *d, const double pos[3]);   /* O(log P) */
```

### Estudiante B — La medición que decide si seguir

**Antes de reescribir la migración**, hay que confirmar que ORB efectivamente
achica el LET. Se puede medir sin tocar nada del pipeline:

1. Correr el particionamiento ORB para obtener las P cajas.
2. Para cada partícula local, calcular su dueño ORB (sin migrar).
3. Construir el árbol local **como si** las partículas ya estuvieran repartidas
   por ORB, y correr `select_nodes` contra las cajas ORB.
4. Contar cuántas entradas se exportarían.

Es un harness de medición, no una implementación: da el número de fantasmas que
tendría el LET con ORB en media jornada, en vez de al final de la semana.

**Criterio de continuación (GO/NO-GO del Día 1):** el volumen estimado a N=100K,
P=4 tiene que bajar de los 61.056 actuales a **menos de 25.000**. Si no baja, el
diagnóstico está mal y hay que parar y replantear, no seguir reescribiendo.

---

## Día 2: Migración por posición

### Estudiante A

El cambio es acotado: en `migrate_particles`, `dest` pasa de búsqueda binaria
sobre splitters a descenso por el árbol ORB.

```c
dest[i] = domain_owner_pos(d, (*p)[i].pos);   /* antes: domain_owner(d, key) */
```

Todo lo demás (`Alltoall` de conteos, `Alltoallv`, gestión de capacidad) queda
igual.

### ⚠ Morton cambia de rol, y hay invariantes que dependían del viejo

Morton deja de determinar la propiedad y queda **solo para localidad de caché**:
ordenar el arreglo local mejora la coherencia en la construcción del árbol y en el
recorrido. Hay que revisar explícitamente lo que asumía lo contrario:

| Lugar | Qué asumía | Qué hacer |
|---|---|---|
| `replicate_particles` (doc) | que el `Allgatherv` sale globalmente ordenado por Morton | **Deja de ser cierto.** Corregir el comentario. Los tests que usan `replicate` emparejan por `id`, así que no dependen del orden — verificar que sea así en los 5 que lo usan. |
| `domain_partition` (doc) | arreglo local ordenado por clave para la búsqueda binaria | La función se retira; ORB no lo necesita. |
| `Domain::gmin/gmax` congelado | claves comparables entre procesos | **Sigue haciendo falta**, ahora por el LET: `octree_build_box` exige la misma caja raíz en todos los procesos. No relajarlo. |
| Test 12 (splitters válidos) | existen splitters | Reescribir (ver Día 3). |

**El bug de `morton_encode` sin acotar (semana 3) sigue siendo relevante:** las
claves se siguen calculando con el box congelado para ordenar localmente.

### Estudiante B — Bandera de comparación

```
--decomp orb|morton      # default: orb
```

Conservar el camino Morton funcionando, exactamente por la misma razón por la que
la semana 4 conservó `--exchange replicate`: **sin poder correr los dos no hay
forma de medir la mejora**, y la comparación es el resultado central del informe.
Fue lo que permitió construir la tabla (a) de la semana 4.

---

## Día 3: Conectar el LET y verificar

### Ambos

`let.c` **no necesita cambios de lógica**. `let_gather_domain_boxes` sigue
publicando el AABB de las partículas locales, que con ORB ya es compacto.

> **Alternativa a evaluar, no a asumir:** el árbol ORB ya le da a cada proceso la
> caja de todos los demás, así que el `Allgather` de cajas podría eliminarse. Pero
> el AABB de las partículas es **más ajustado** que la caja ORB (que se extiende
> hasta los bordes del dominio global, inflado por los outliers de Plummer).
> Medir las dos; probablemente gane el AABB, y el Allgather de 6 doubles no es un
> costo que valga la pena pelear.

### Tests

**Test 12 — reescrito para ORB.** Ya no hay splitters:
- las cajas de los P procesos son **disjuntas** (intersección de volumen nula);
- **cubren** la caja global;
- toda partícula local cae dentro de la caja de su proceso;
- desbalance de trabajo < 1.05.

**Test 14 (θ=0 bit a bit) debe seguir pasando.** Es el guardián de regresión más
fuerte que hay: si la migración por posición pierde o duplica una partícula, o si
la caja raíz deja de ser común, el test 14 lo detecta con error distinto de cero.
Correrlo apenas compile la migración nueva, antes que ningún benchmark.

**Test 16 pasa de ser una falla documentada a ser el criterio de aceptación:**
exponente k < 0.9. Es el test que define si la semana sirvió.

**Test 19 (nuevo) — Volumen del LET vs P.** A N fijo, medir `n_ghost` con
P=2,4,8,16. Con dominios compactos debe crecer como superficie, no como volumen.
Criterio: `n_ghost/n_local` **no crece** proporcionalmente a P (con Morton se medió
1,00 → 2,44 → 4,36 para P=2,4,8, que es aproximadamente lineal en P).

### Punto de control del Día 3

Si el test 16 no da k < 0.9 acá, **quedan dos días y hay que decidir**: o se
diagnostica con la misma disciplina de la semana 4 (medir, no suponer), o se
congela y se escribe el informe con dos resultados negativos bien documentados.
No entrar al Día 4 sin haber tomado la decisión.

---

## Día 4: Experimentos

Todo en `pcunix*` (`module load mpi/mpich-x86_64`, `export FI_PROVIDER=tcp`,
`OMP_PROC_BIND=close`, `OMP_PLACES=cores`).

> **Deuda acumulada de dos semanas.** Las tablas de las semanas 3 y 4 salen de un
> portátil con 4 núcleos performance + 4 efficiency. La tabla híbrida de la semana
> 4 salió tan ruidosa que quedó marcada como inservible. **Re-medir las tres
> semanas en pcunix no es opcional**, y es lo primero del Día 4, no lo último.

| # | Experimento | Qué debe mostrar |
|---|---|---|
| a | **n_ghost vs N**, ORB vs Morton, P=4 | el exponente k baja de 0,92 a < 0,9. **El gráfico de la semana.** |
| b | **n_ghost vs P**, N fijo | crecimiento de superficie, no de volumen |
| c | **Desglose por fase**, ORB vs Morton vs replicación, P=1..16 | `tree_time` deja de crecer con P |
| d | **Escalabilidad fuerte**, N=100K y 500K | speedup > 3× con P=4 |
| e | **Escalabilidad débil**, N/P=25K | eficiencia > 60% con P=8 (semana 4: 5,4%) |
| f | **Híbrido**, P×T constante | rehacer la tabla que quedó inservible |
| g | **Plummer vs uniform**, ORB | la brecha 61.056 / 22.029 debería cerrarse |
| h | **Distribuido**, `-hosts pcunix40,pcunix42` | deuda desde la semana 3 |

El experimento **(g)** es la verificación más directa del diagnóstico: si la causa
era la geometría de los dominios sobre distribuciones concentradas, con ORB la
diferencia entre Plummer y uniform tiene que achicarse mucho.

---

## Día 5: Informe final

`docs/final-report.md`. Estructura: problema y algoritmo → arquitectura híbrida →
descomposición del dominio (Morton → ORB, con el porqué) → LET → balance dinámico
→ validación → resultados → análisis crítico → conclusiones.

El eje narrativo del informe ya está escrito por las mediciones, y es lo más
valioso que tiene el trabajo:

1. **Semana 3:** se esperaba que dominara la comunicación. Dominaba la
   construcción redundante del árbol (`Allgatherv` = 1,3%).
2. **Semana 4:** se esperaba que el LET redujera el árbol local a N/P. Lo redujo a
   ~N/2, porque los dominios Morton no son compactos.
3. **Semana 5:** ORB.

Las tres veces la hipótesis razonable estaba equivocada y la instrumentación lo
mostró. Vale la pena decir explícitamente que **el LET de la semana 4 era
correcto** —validado bit a bit con θ=0— **y aun así no servía para lo que se lo
puso**: correcto y efectivo son cosas distintas, y solo la medición las distingue.

---

## Criterios de aceptación

| Criterio | Umbral | Línea base (semana 4) |
|---|---|---|
| **`n_ghost` ~ N^k** | **k < 0.9** | **k = 0.92** |
| **`tree_time` P=8 vs P=1** | **no crece** | **4,6×** |
| `n_ghost/n_local` a N=100K, P=4 | < 1.0 | 2,44 |
| LET con θ=0 vs BH secuencial | bit a bit | OK (mantener) |
| Cajas ORB disjuntas y cubrientes | exacto | — (test 12 nuevo) |
| Checksum Σid, Σid² tras migración por posición | exacto | OK (mantener) |
| Desbalance de trabajo | < 1.05 | 1.0013 (mantener) |
| Energía, 200 pasos | drift < 1% | 1.21e-04 (mantener) |
| Speedup MPI, P=4 | > 3× | no alcanzado |
| Eficiencia débil, P=8 | > 60% | 5,4% |
| Suite completa | **19/19** | 17/18 |

## Riesgos

| Riesgo | Impacto | Mitigación |
|---|---|---|
| **ORB no arregla el problema** | **Alto** | Harness de medición del Día 1 **antes** de reescribir nada; GO/NO-GO explícito |
| Partículas perdidas en la migración por posición | Alto | Test 14 (bit a bit) y test 11 (checksum) apenas compile |
| Discrepancia `<` vs `>=` en el plano de corte | Alto | Mismo predicado en el conteo y en `domain_owner_pos`; test 12 lo cubre |
| Bisección lineal cara | Medio | Medirla en `metrics`; histograma de bins si aparece |
| P no potencia de 2 | Medio | Partición proporcional desde el principio, no como parche |
| Se rompe algo que asumía propiedad por clave | Medio | Tabla de invariantes del Día 2, revisada explícitamente |
| **No queda tiempo para el informe** | **Alto** | Congelar código el Día 3; los Días 4-5 son medir y escribir |
| Las mediciones de pcunix no llegan | Alto | Día 4 empieza por ahí |

## Prioridades si falta tiempo

1. ORB + migración por posición + test 14 en verde (sin esto no hay semana).
2. Test 16 con k < 0.9: el criterio que define el éxito.
3. Experimento (a): n_ghost vs N, ORB vs Morton.
4. Re-medición en pcunix de las semanas 3-5.
5. Escalabilidad fuerte y débil.
6. Experimentos f, g, h.

**El informe final no es sacrificable.** Si el Día 3 ORB no funciona, congelar y
escribir: dos resultados negativos bien medidos y bien explicados valen más que
una tercera reescritura a medio validar.

## Entregables

1. `src/domain.h/.c` reescrito con ORB; `domain_owner_pos`.
2. `migrate_particles` por posición; invariantes de Morton revisadas y documentadas.
3. `--decomp orb|morton` con ambos caminos funcionando.
4. Test 12 reescrito, test 16 como criterio de aceptación, test 19 nuevo → 19/19.
5. `docs/specs/orb.md`; `let.md` y `domain.md` actualizados.
6. `docs/week5-report.md` con los experimentos a-h.
7. `docs/final-report.md`.
