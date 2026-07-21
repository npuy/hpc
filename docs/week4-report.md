# Semana 4 — Reporte: LET, Balance Dinámico y Resultado Negativo

## Resumen

Se implementaron las dos etapas previstas —LET (etapa 8) y rebalanceo dinámico
(etapa 9)— y ambas funcionan y están validadas. El rebalanceo por costo cumple su
criterio con holgura. **El LET es correcto pero no comprime**: el volumen de datos
importados crece linealmente con N, la misma asintótica que la replicación a la
que venía a reemplazar, y el criterio central de la semana no se cumple.

La causa está medida y es un supuesto equivocado del plan de la semana, no un bug:
**un tramo contiguo de claves Morton no es una región compacta del espacio**. La
sección "Por qué el LET no comprime" documenta la evidencia.

## Archivos implementados

```
src/
  let.h / .c          -- NUEVO: seleccion, empaquetado e importacion del LET
  octree.h / .c       -- octree_build_box (caja raiz explicita)
                         octree_insert_particles (insercion incremental)
  barnes_hut.c        -- contador de interacciones -> Particle::work
  particle.h          -- campo work, convencion id == -1 para fantasmas
  mpi_types.c         -- 4o bloque en el tipo derivado (work)
  domain.h / .c       -- domain_partition_ex (biseccion ponderada por costo)
  metrics.h / .c      -- fases let_time y ghost_sum
  main.c              -- build_step_tree, flags --exchange/--balance/--theta-let
  validation.h / .c   -- tests 14-18
docs/specs/
  let.md, balance.md
```

## Validación: 17/18

| Test | Descripción | Criterio | Resultado | Valor |
|---|---|---|---|---|
| 1-13 | Suite de semanas 1-3 | — | **PASS** | sin regresión |
| 14 | LET con θ=0 vs BH secuencial | bit a bit | **PASS** | **0.0 exacto** |
| 15 | Precisión del LET a θ=0.5 | < 1e-3 y < 1.5× err. BH | **PASS** | ver abajo |
| 16 | Volumen del LET (crecimiento) | n_ghost ~ N^k, k < 0.9 | **FAIL** | **k = 0.92** |
| 17 | Rebalanceo por costo | desbalance < 1.05 | **PASS** | 1.2072 → **1.0013** |
| 18 | Energía con LET, 200 pasos | drift < 1% | **PASS** | 1.21e-04 |

El **test 14** es el que garantiza que la selección no se saltea ninguna
interacción: con θ=0 el criterio de exportación no acepta ningún resumen, el LET
degenera en réplica completa, el árbol fusionado tiene la misma estructura y el
mismo orden de sumas que el secuencial, y la igualdad vuelve a ser bit a bit. Es
el reemplazo del test 13, cuyo criterio (error exactamente 0.0) deja de tener
sentido con θ>0.

El **test 16 falla a propósito y no hay que "arreglarlo" bajando el umbral.** Es
la medición del resultado negativo de la semana.

### Sobre el test 15

```
LET vs Barnes-Hut secuencial:              0.000000e+00   (umbral 1e-3)
Barnes-Hut secuencial vs directo O(N^2):   3.081068e-02
LET vs directo O(N^2):                     3.081068e-02   (umbral 4.62e-02)
```

El error del LET contra el secuencial da **exactamente cero**, lo que a primera
vista parece un resultado excelente. No lo es: es el primer síntoma del problema.
A N=2000 con P=4 el LET casi no resume nada, así que cada proceso termina viendo
prácticamente todas las partículas y el cálculo es idéntico al replicado. El error
de reagrupamiento que el plan anticipaba (un fantasma reinsertado puede quedar
agrupado bajo un ancestro y aproximarse dos veces) **no llegó a manifestarse
porque casi no hay fantasmas resumidos**.

## Por qué el LET no comprime

### La medición

Fantasmas importados por proceso, P=4, Plummer:

| N | n_local | fantasmas | ratio |
|---|---|---|---|
| 12.500 | 3.125 | 8.324 | 2,66 |
| 25.000 | 6.250 | 16.271 | 2,60 |
| 50.000 | 12.500 | 33.453 | 2,67 |
| 100.000 | 25.000 | 61.056 | 2,44 |
| 200.000 | 50.000 | 114.289 | 2,28 |
| 400.000 | 100.000 | 243.456 | 2,43 |

**El ratio es constante**, o sea que n_ghost ∝ N. El exponente medido por el test
16 es k=0.92, indistinguible de lineal. Un LET que funciona debería mostrar un
ratio que cae al crecer N, porque el término de superficie pierde peso frente al
de volumen.

### La pista que orientó el diagnóstico

`--theta-let` no tiene **ningún** efecto sobre el volumen exportado:

| θ_export | 0,5 | 1,0 | 2,0 | 4,0 |
|---|---|---|---|---|
| fantasmas | 61.056 | 60.670 | 66.438 | 63.195 |

Si el criterio de apertura estuviera gobernando la selección, multiplicar θ por 8
tendría que desplomar el volumen. Que no cambie nada significa que **la rama que
acepta el resumen casi nunca se toma**: se desciende siempre hasta las hojas y se
exportan partículas reales.

### La causa

Instrumentando las cajas de dominio que publican los procesos (N=100K, P=4,
lado del dominio global L=410,6):

```
caja[0] lados=(113.0 103.0 151.9)   vol rel=0.0255   envio=0
caja[1] lados=(  3.2   3.2   4.8)   vol rel=0.0000   envio=15045
caja[2] lados=( 86.2  79.9  87.4)   vol rel=0.0087   envio=24998
caja[3] lados=(402.5 293.2 261.3)   vol rel=0.4456   envio=24998
```

La caja del proceso 3 cubre el **44,6% del volumen del dominio**. Contra un
destino así, la distancia mínima `d_min` de casi cualquier nodo a la caja es cero,
el MAC nunca se satisface y se exporta todo sin resumir: `envio=24998` son
exactamente las 25.000 partículas locales del proceso 0.

El supuesto que falla viene del plan de la semana: **se dio por sentado que un
tramo Morton contiguo describe una región compacta, y no es así.** Dos efectos se
suman, y el segundo es el grave:

1. **El AABB no es robusto a la cola pesada de Plummer.** El core mide ~3 unidades
   y el bounding box global mide 410, estirado por un puñado de partículas a radio
   grande. Basta una para inflar la caja de un proceso.

2. **El pico de densidad cae justo en la esquina donde se tocan los octantes de
   nivel 1.** El box global está centrado en el centro de masa, que es donde
   Plummer concentra casi toda su masa. Ahí se juntan los 8 octantes de primer
   nivel, así que las claves Morton del core se reparten por *todo* el espacio de
   claves. Cualquier corte por percentil de conteo termina dando a cada proceso
   una mezcla de core y halo, y los dominios se interpenetran geométricamente.

El segundo punto es el que hace que no haya arreglo barato. Si el dominio del
proceso 3 contiene genuinamente partículas pegadas al core, el detalle del core
hay que mandárselo: no es que el descriptor sea pesimista, es que el dominio está
efectivamente esparcido. Cambiar el AABB por un conjunto de celdas describiría
mejor la región pero no cambiaría la conclusión.

**Evidencia de apoyo:** con `--init uniform`, donde no hay pico central, el mismo
caso baja de 61.056 a **22.029** fantasmas. La patología es de la distribución
centralmente concentrada combinada con descomposición por rangos Morton, no del
código.

### Qué lo arreglaría

Dominios compactos por construcción, es decir **ORB (bisección recursiva
ortogonal)**: cortar el espacio por planos balanceando trabajo, de modo que cada
proceso reciba una caja. Con dominios que son cajas, el criterio del AABB pasa a
ser exacto y ajustado, y el LET recupera su comportamiento esperado. La maquinaria
de bisección por histograma de `domain.c` se reusa casi tal cual (bisectando una
coordenada en vez del espacio de claves), pero cambia `domain_owner`, la migración
y el test 12. Queda como el trabajo futuro prioritario.

## Resultados de rendimiento

> **Entorno:** mismo portátil Apple Silicon de la semana 3 (4 núcleos performance +
> 4 efficiency, Open MPI 5.0.9, clang + libomp). **No** son las máquinas `pcunix`.
> Con P=8 se sobresuscribe y se usan los núcleos E, que rinden ~21% de un núcleo P,
> así que **las filas de P=8 no son confiables** y las de configuración híbrida son
> directamente ruidosas. Todo esto hay que re-medirlo en pcunix.

### a) Replicación vs LET (N=100K, 10 pasos, 1 hilo)

| Config | Total (s) | Árbol (s) | Fuerzas (s) | Réplica (s) | LET (s) | Migración (s) |
|---|---|---|---|---|---|---|
| P=1 replicate | 3.618 | 0.247 | 3.032 | 0.002 | — | 0.052 |
| P=1 let | 3.612 | 0.256 | 2.957 | — | 0.003 | 0.052 |
| P=2 replicate | 3.847 | 0.311 | 3.035 | 0.008 | — | 0.146 |
| P=2 let | 3.894 | 0.362 | 2.999 | — | 0.064 | 0.148 |
| P=4 replicate | 4.725 | 0.762 | 3.399 | 0.024 | — | 0.155 |
| P=4 let | **4.188** | 0.736 | 2.903 | — | 0.118 | 0.090 |
| P=8 replicate | 5.392 | 1.947 | 2.672 | 0.110 | — | 0.248 |
| P=8 let | **4.395** | **1.184** | 2.151 | — | 0.297 | 0.304 |

El LET **sí mejora**, pero parcialmente: 18% menos tiempo total y 39% menos tiempo
de árbol con P=8. Lo que no logra es el objetivo de la semana:

| | P=1 → P=8 |
|---|---|
| `tree_time` con replicación | 0.247 → 1.947 s (**7,9×**) |
| `tree_time` con LET | 0.256 → 1.184 s (**4,6×**) |

**El criterio era que dejara de crecer.** Sigue creciendo, porque el árbol se
construye sobre `n_local + n_ghost` y n_ghost ≈ 2,5·n_local: el árbol local tiene
del orden de la mitad de las partículas globales en vez de N/P.

### b) Escalabilidad débil (N/P = 25K, 1 hilo)

| P | N | replicate (s) | let (s) | árbol repl. (s) | árbol let (s) |
|---|---|---|---|---|---|
| 1 | 25.000 | 0.734 | 0.671 | 0.053 | 0.059 |
| 2 | 50.000 | 1.637 | 1.650 | 0.124 | 0.127 |
| 4 | 100.000 | 4.334 | 4.337 | 0.732 | 0.669 |
| 8 | 200.000 | 15.054 | **12.530** | 5.276 | **3.267** |

Eficiencia débil a P=8: 4,9% con replicación, **5,4% con LET**. La mejora es real
pero marginal, y ninguna de las dos escala. Es coherente con n_ghost ∝ N: al
mantener N/P fijo el trabajo por proceso sigue creciendo con P.

### c) Rebalanceo por costo (N=100K, LET, 20 pasos)

| P | Reparto | Total (s) | Desbal. partículas | Desbal. **trabajo** |
|---|---|---|---|---|
| 4 | count | 8.948 | 1.0015 | 1.1382 |
| 4 | **work** | 8.970 | 1.0723 | **1.0473** |
| 8 | count | 9.815 | 1.0030 | 1.2469 |
| 8 | **work** | **9.142** | 1.2142 | **1.1245** |

**Este es el resultado limpio de la semana.** El reparto por costo baja el
desbalance de trabajo en los dos casos, y a P=8 se traduce en 7% menos tiempo
total. Y se observa exactamente la inversión que el plan anticipaba como señal de
funcionamiento: el desbalance de **partículas empeora** (1.003 → 1.214) mientras
el de **trabajo mejora** (1.247 → 1.124), porque los procesos que reciben zonas
densas se quedan con menos partículas.

En el test 17, con reparticiones más frecuentes, el efecto es más nítido todavía:
**1.2072 → 1.0013**.

### d) Precisión vs θ (N=20K, P=4, LET)

| θ | Tiempo (s) | Energy drift |
|---|---|---|
| 0,2 | 7.202 | 3.55e-06 |
| 0,5 | 0.744 | 2.60e-05 |
| 0,8 | 0.291 | 2.34e-05 |
| 1,0 | 0.299 | 2.80e-05 |

El compromiso costo/precisión se comporta como corresponde a Barnes-Hut y la
energía se conserva bien en todo el rango.

### e) Configuración híbrida (P×T=8, N=100K, LET)

| Config | Total (s) | Árbol (s) | Fuerzas (s) |
|---|---|---|---|
| 1×8 | 4.413 | 0.300 | 3.615 |
| 2×4 | 5.122 | 0.379 | 4.277 |
| 4×2 | **4.327** | 0.710 | 2.951 |
| 8×1 | 5.351 | 1.087 | 3.003 |

Los números no son monótonos (2×4 sale peor que 1×8 y que 4×2), lo que delata
ruido de medición en un portátil con núcleos heterogéneos. **No se puede concluir
nada de esta tabla y hay que rehacerla en pcunix.** La pregunta que debía
responder —si el LET invierte la conclusión de la semana 3 de que conviene
minimizar procesos— queda abierta.

## Criterios de aceptación

| Criterio | Umbral | Estado | Valor |
|---|---|---|---|
| LET con θ=0 == BH secuencial | bit a bit | **OK** | 0.0 exacto |
| LET (θ=0.5) vs replicado | < 1e-3 | **OK** | 0.0 |
| Error LET / error BH secuencial | < 1.5× | **OK** | 1.00× |
| Σ masa (locales + fantasmas) | < 1e-12 | **OK** | 4.9e-14 |
| Checksum Σid, Σid² (sin fantasmas) | exacto | **OK** | test 11 |
| Energía con LET, 200 pasos | < 1% | **OK** | 1.21e-04 |
| **`tree_time` P=8 vs P=1** | **no crece** | **NO** | **4,6×** (repl.: 7,9×) |
| Speedup MPI, P=4 | > 3× | **NO** | no medible aquí |
| Eficiencia débil, P=8 | > 60% | **NO** | 5,4% |
| Desbalance de trabajo, `--balance work` | < 1.05 | **OK** | 1.0013 (test 17) |
| **n_ghost, crecimiento con N** | **sublineal** | **NO** | **k = 0.92** |
| Sin fugas ni deadlocks | limpio | **OK** | solo colectivos |
| Suite de validación | 18/18 | **NO** | **17/18** |

**9 de 13 criterios se cumplen.** Los cuatro que no comparten una única causa: el
LET no comprime, así que el árbol local sigue siendo grande y todo lo que dependía
de que se achicara queda sin cumplirse.

## Estado del código

Ambos esquemas conviven y se eligen por CLI, que es lo que permitió la comparación
directa de la sección (a):

```bash
--exchange let|replicate   # LET (default) o replica global de la semana 3
--balance work|count       # reparto por costo (default) o por conteo
--theta-let X              # theta del criterio de exportacion (default: = --theta)
```

`replicate_particles` se conserva a propósito: es el oráculo de los tests y la
línea base de los gráficos.

## Pendiente

- **ORB para dominios compactos.** Es el trabajo prioritario y lo único que
  arreglaría el resultado negativo de fondo.
- **Re-medir todo en pcunix**, incluida la semana 3. Las tablas de este reporte y
  las de la semana 3 salen de un portátil con núcleos heterogéneos; sin línea base
  homogénea las comparaciones no se sostienen. La tabla (e) es directamente
  inservible.
- **Corrida distribuida en varias máquinas** (`-hosts`): sigue sin verificarse,
  arrastrada desde la semana 3.
- **Informe final** (`docs/final-report.md`).

## Lectura de la semana

Es la segunda semana seguida en que la medición contradice a la hipótesis del
plan, y en las dos el mecanismo fue el mismo: una suposición razonable sobre dónde
estaba el costo, que nadie había verificado.

- **Semana 3:** se esperaba que dominara la comunicación. Dominaba la construcción
  redundante del árbol (el `Allgatherv` era el 1,3%).
- **Semana 4:** se esperaba que el LET redujera el árbol local a N/P. Lo reduce a
  ~N/2, porque los dominios Morton no son compactos sobre una distribución
  centralmente concentrada.

Lo que hizo visible el problema no fue leer el código sino dos mediciones baratas:
que `--theta-let` no cambiara nada, y el volumen de las cajas de dominio. Vale la
pena señalarlo en el informe final: el LET *funciona* —está validado bit a bit con
θ=0— y aun así no sirve para lo que se lo puso. Correcto y efectivo son cosas
distintas, y solo la instrumentación distingue una de la otra.
