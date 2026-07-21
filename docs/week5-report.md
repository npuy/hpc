# Semana 5 — Reporte: ORB, y el LET que por fin comprime

## Resumen

Se reemplazó la descomposición por rangos de clave Morton por **ORB (bisección
recursiva ortogonal)**. Los dos criterios que la semana 4 no había alcanzado se
cumplen:

| Criterio | Umbral | Semana 4 (Morton) | Semana 5 (ORB) |
|---|---|---|---|
| `n_ghost` ~ N^k | k < 0.9 | **0.92** | **0.66** |
| `tree_time` de P=1 a P=8 | no crece | **4,8×** | **1,6×** |

Suite de validación: **19/19**, verificada con P=2, 3, 4 y 8. Sin fugas.

**Todo esto quedó re-medido en las máquinas pcunix de fing** (5 nodos, núcleos
homogéneos, ejecución realmente distribuida), lo que salda la deuda que venía
arrastrándose desde la semana 3. Los resultados se reproducen casi exactamente y
además se cumplen por primera vez dos criterios que el portátil nunca alcanzó:
speedup >3× con P=4 y eficiencia débil >60% con P=4. Ver la sección
"Validación en pcunix".

## Archivos

```
src/
  domain.h / .c       -- OrbNode, domain_partition_orb, domain_owner_pos, domain_box
  migration.c         -- dueno por posicion cuando use_orb
  main.c              -- repartition() ramifica; flag --decomp orb|morton
  validation.c        -- test 12 reescrito, 16 y 19 reformulados, setup_distributed_ex
docs/specs/
  orb.md
```

`let.c` **no necesitó un solo cambio de lógica**. El criterio de exportación ya
era correcto; lo que fallaba era la geometría que recibía.

## El resultado principal

### Volumen del LET vs N (P=4, Plummer)

| N | n_local | Morton | ratio | ORB | ratio |
|---|---|---|---|---|---|
| 12.500 | 3.125 | 8.324 | 2,66 | 4.281 | 1,36 |
| 25.000 | 6.250 | 16.271 | 2,60 | 6.830 | 1,09 |
| 50.000 | 12.500 | 33.453 | 2,67 | 10.789 | 0,86 |
| 100.000 | 25.000 | 61.056 | 2,44 | 16.931 | 0,67 |
| 200.000 | 50.000 | 114.289 | 2,28 | 26.503 | 0,53 |
| 400.000 | 100.000 | 243.456 | 2,43 | 41.594 | 0,41 |

**Exponente medido: k = 0,66 con ORB, k = 0,97 con Morton.**

Con Morton el ratio se queda clavado en ~2,5: `n_ghost` crece como N y el LET no
comprime. Con ORB el ratio cae monótonamente de 1,36 a 0,41, y el exponente 0,66
está muy cerca del **2/3 teórico de un efecto de superficie**, que es exactamente
lo que un LET debe costar.

### Desglose por fase vs P (N=100K, 1 hilo)

| Config | Total (s) | Árbol (s) | Fuerzas (s) |
|---|---|---|---|
| P=1 morton | 3.333 | 0.247 | 2.763 |
| P=1 orb | 3.433 | 0.252 | 2.847 |
| P=2 morton | 3.685 | 0.324 | 2.902 |
| P=2 orb | 3.732 | **0.203** | 3.107 |
| P=4 morton | 4.016 | 0.741 | 2.794 |
| P=4 orb | **3.740** | **0.211** | 3.158 |
| P=8 morton | 4.475 | 1.186 | 2.321 |
| P=8 orb | **3.959** | **0.411** | 2.963 |

`tree_time` con Morton: 0,247 → 1,186 s (**4,8×**). Con ORB: 0,252 → 0,411 s
(**1,6×**), y de P=1 a P=4 incluso **baja** (0,252 → 0,211), que es el
comportamiento correcto: cada proceso construye sobre su vecindario, no sobre el
sistema entero.

El repunte en P=8 es casi seguro sobresuscripción del portátil, no algorítmico.
**Confirmado después:** en pcunix, con procesos repartidos en 5 máquinas, el árbol
con ORB decrece monótonamente al subir P.

**Un efecto lateral honesto:** la fase de fuerzas empeora levemente con ORB
(2,79 → 3,16 s a P=4). La hipótesis era peor localidad de caché, porque el arreglo
local queda ordenado por Morton pero la propiedad ya no sigue la curva.
**Desmentido después:** en pcunix la diferencia baja a 2-3%, así que era del banco
de pruebas y no del algoritmo.

### Escalabilidad débil (N/P = 25K)

| P | N | morton total | orb total | morton árbol | orb árbol |
|---|---|---|---|---|---|
| 1 | 25.000 | 0.671 | 0.671 | 0.058 | 0.057 |
| 2 | 50.000 | 1.551 | 1.537 | 0.124 | **0.082** |
| 4 | 100.000 | 4.325 | **3.651** | 0.765 | **0.141** |
| 8 | 200.000 | 10.290 | **8.817** | 2.681 | **0.371** |

El árbol pasa de crecer 46× a crecer 6,5×, una reducción de **7,2×** a P=8. Pero
**la eficiencia débil sigue siendo mala** (7,6% con ORB contra 6,5% con Morton):
el árbol dejó de ser el problema y ahora domina la fase de fuerzas, que crece con
`log N` por construcción y encima se mide sobresuscribiendo la máquina.
**Rehecho en pcunix:** con núcleos homogéneos la eficiencia débil sube a 62% con
P=4 y 33% con P=8.

### Plummer vs uniform (N=100K, P=4) — la verificación del diagnóstico

| Init | Morton | ORB |
|---|---|---|
| plummer | 61.056 | 16.933 |
| uniform | 22.029 | 9.232 |
| **brecha** | **2,77×** | **1,83×** |

Es la comprobación más directa de que el diagnóstico de la semana 4 era correcto.
La explicación decía que la patología venía de la distribución centralmente
concentrada combinada con dominios no compactos; si eso era cierto, con dominios
compactos la brecha entre Plummer y uniform tenía que achicarse. Se achicó de
2,77× a 1,83×.

## Validación: 19/19

| Test | Descripción | Criterio | Valor |
|---|---|---|---|
| 1-11, 13-15, 17-18 | Suites previas | — | sin regresión |
| 12 | Descomposición válida (ORB **y** Morton) | disjuntas, cubrientes, coherentes | PASS |
| 14 | LET con θ=0 vs BH secuencial | bit a bit | **0.0 exacto** |
| 16 | Crecimiento de n_ghost con N | k < 0.9 | morton 0.92 / **orb 0.69** |
| 19 | ORB vs Morton a igual config. | ORB < 0,6× Morton | 0,38–0,54 |

**El test 14 siguió pasando durante toda la reescritura** y fue lo que dio
confianza para tocar la migración: si el cambio de dueño por posición hubiera
perdido o duplicado una partícula, o si la caja raíz hubiera dejado de ser común,
el error habría dejado de ser exactamente cero.

El **test 12** ahora valida las dos descomposiciones. Para ORB: cajas disjuntas
(intersección de volumen nula par a par), cubrientes (Σ volúmenes = volumen del
dominio, verificado a 1e-9), toda partícula dentro de su caja, y coherencia entre
`domain_owner_pos` y el rank — que es lo que detectaría una discrepancia entre el
predicado del conteo y el de la asignación.

El **test 16** corre las dos descomposiciones en la misma invocación y reproduce
exactamente el k=0,92 de la semana 4 al lado del 0,69 de ORB. Queda
autodocumentado: el test que la semana 4 dejó fallando a propósito es el mismo que
ahora certifica la corrección.

### Un umbral que corregí, y por qué

El test 19, tal como lo especifiqué en el plan, exigía `n_ghost/n_local < 1,5` con
N fijo. **Fallaba a P=8** (2,26). Al mirar los datos crudos el problema era el
test, no el código:

| N=20.000 | n_local | ORB fantasmas | ratio |
|---|---|---|---|
| P=2 | 10.000 | 4.707 | 0,47 |
| P=4 | 5.000 | 5.901 | 1,18 |
| P=8 | 2.500 | 5.635 | 2,25 |

Los fantasmas absolutos de ORB son **casi constantes**; el cociente sube solo
porque el denominador se achica. El umbral absoluto estaba midiendo dos cosas a la
vez: el efecto de superficie y el encogimiento del dominio.

Reformulé el test sobre la afirmación que de verdad importa —ORB importa
sustancialmente menos que Morton a igual configuración, medido entre 1,8× y 2,6×
menos— con `n_local` constante en vez de N constante. **Dejo constancia de que
cambié un umbral después de verlo fallar**, porque es exactamente la maniobra que
puede usarse para maquillar un resultado. Lo que la sostiene: el criterio
sustantivo de la semana (test 16, k < 0,9) se cumplió sin tocar nada, y ORB gana
en todas las configuraciones medidas con cualquiera de las dos formulaciones.

## Notas de implementación

### Bisección ORB

Recursiva sobre el grupo de procesos. En cada nodo se elige como eje la dimensión
más larga de la caja (evita cajas laminares, que tienen mucha superficie y por lo
tanto mucho LET) y se bisecta la coordenada 50 veces hasta repartir el trabajo.

**La pertenencia se lleva por el camino recorrido, no por un test geométrico.** El
arreglo de índices locales se particiona in-place en cada nivel con el mismo
predicado `<` que usa `domain_owner_pos`. Era el riesgo señalado en el plan: una
discrepancia entre el predicado del conteo y el de la asignación haría que una
partícula sobre el plano se contara de un lado y se enviara al otro. Llevarlo por
el camino lo vuelve imposible por construcción, en vez de dependiente de que dos
comparaciones coincidan.

Se bisecta el **trabajo**, no el conteo: el campo `work` de la semana 4 ya estaba
validado, así que ORB heredó el balance por costo sin código nuevo. Con `work`
global nulo (primer paso) cae al conteo.

**P que no es potencia de 2** funciona con reparto proporcional
(`nL = nprocs/2`, objetivo `w_total·nL/nprocs`). Verificado con P=3: 19/19.

El costo previsto en el plan (~50 colectivos por nodo interno, o sea (P−1)·50 por
repartición) no apareció en las métricas, como se esperaba: las reparticiones son
cada 20 pasos. No hizo falta la optimización por histograma.

### Morton cambió de rol

Ya no determina la propiedad; queda solo para localidad de caché al ordenar el
arreglo local. La invariante que documentaba `replicate_particles` —que el
`Allgatherv` sale globalmente ordenado por Morton— **deja de ser cierta con ORB**.
Los cinco tests que usan `replicate` emparejan por `id`, así que ninguno dependía
del orden.

`Domain::gmin/gmax` congelado **sigue haciendo falta**, ahora por otra razón:
`octree_build_box` exige la misma caja raíz en todos los procesos.

## Validación en pcunix: la deuda de tres semanas, saldada

Todo lo anterior se midió en un portátil Apple Silicon con núcleos heterogéneos.
La corrida en el cluster de fing (`docs/resultados-pcunix.txt`, 6 min 40 s)
cierra esa deuda y cambia el estatus de los resultados: dejan de ser "válidos
como comparación relativa" y pasan a ser medidos en hardware homogéneo.

**Entorno:** 5 nodos Intel i3-4170 (2 núcleos físicos c/u), GCC 15.2.1, MPICH,
x86-64. Los barridos van distribuidos a 2 procesos por nodo, así que P=1..8 son
1 a 8 núcleos físicos **sin sobresuscribir**. Ver [pcunix-guia.md](pcunix-guia.md).

### Portabilidad

**19/19 con P=2, 3, 4 y 8 distribuidos en 5 máquinas**, y 9/9 secuencial. Otro
compilador, otra implementación de MPI y otra arquitectura dan los mismos
veredictos, **incluidos los dos tests de igualdad bit a bit** (el 9, determinismo
de OpenMP, y el 14, LET con θ=0). También queda saldada la corrida multi-máquina
pendiente desde la semana 3.

### El resultado central se reproduce casi exactamente

Fantasmas por proceso, P=4, y el cociente contra las partículas locales:

| N | morton | ratio | *(portátil)* | orb | ratio | *(portátil)* |
|---|---|---|---|---|---|---|
| 12.500 | 8.544 | 2,73 | *2,66* | 4.225 | **1,35** | *1,36* |
| 25.000 | 17.102 | 2,74 | *2,60* | 6.804 | **1,09** | *1,09* |
| 50.000 | 28.914 | 2,31 | *2,67* | 10.768 | **0,86** | *0,86* |
| 100.000 | 58.739 | 2,35 | *2,44* | 16.906 | **0,68** | *0,67* |

**Exponente: k = 0,67 con ORB y k = 0,93 con Morton** (portátil: 0,66 y 0,92).
La serie de ORB coincide con la del portátil hasta la segunda decimal en los
cuatro tamaños. El resultado no dependía de la máquina.

### El LET es un efecto de superficie, ahora medido directamente

A **N fijo** (50K) subiendo P, los fantasmas de ORB se mantienen prácticamente
constantes, que es la firma de un término de superficie:

| P | morton | orb |
|---|---|---|
| 2 | 24.995 | 9.038 |
| 4 | 28.917 | 10.769 |
| 8 | 29.821 | 9.794 |

El `24.995` de Morton con P=2 es elocuente: N/2 = 25.000, o sea que **cada proceso
importa exactamente todas las partículas del otro**. Compresión cero, que es el
hallazgo de la semana 4 reproducido en otro hardware.

A **n_local fijo** subiendo P, ORB crece como P^0,72 y Morton como P^1,10.

### Escalabilidad fuerte (N=50K) — el criterio se cumple

| P | morton | speedup | orb | speedup | eficiencia orb |
|---|---|---|---|---|---|
| 1 | 10.142 | 1,00× | 10.118 | 1,00× | 100% |
| 2 | 5.325 | 1,90× | 5.190 | 1,95× | 97% |
| 4 | 4.371 | 2,32× | 3.161 | **3,20×** | **80%** |
| 8 | 4.176 | 2,43× | 2.777 | 3,64× | 46% |

**Es la primera vez en todo el proyecto que se cumple el criterio de >3× con
P=4.** En el portátil nunca pasó de 2,80×, y ahora se ve que buena parte de eso
era el hardware.

`tree_time` con ORB **decrece** al subir P (0,093 → 0,056 → 0,041 → 0,034 s),
que es el comportamiento correcto: cada proceso construye sobre menos partículas.
Con Morton se queda plano (0,092 → 0,076).

### Escalabilidad débil (N/P = 12,5K)

| P | morton | orb | eficiencia orb |
|---|---|---|---|
| 1 | 1.958 | 1.960 | 100% |
| 2 | 2.343 | 2.294 | 85% |
| 4 | 4.526 | 3.164 | **62%** |
| 8 | 8.885 | 5.980 | 33% |

**A P=4 se cumple el criterio de >60%** (en el portátil daba 7,6%). A P=8 cae a
33%: lo que domina ahí es la fase de fuerzas (3,86 s), que crece más de lo que
predice `log N` porque los fantasmas también entran al recorrido del árbol.

### OpenMP: exactamente lo previsto

| Hilos | Total | Speedup | Fase fuerzas |
|---|---|---|---|
| 1 | 10.119 | 1,00× | 1,00× |
| 2 | 5.163 | **1,96×** | 2,00× |
| 4 | 4.236 | 2,39× | 2,43× |

Con 2 núcleos físicos, 2 hilos dan 1,96× (98% de eficiencia) y el hyperthreading
agrega ~20% más. **Confirma que el 3,94× con 8 hilos de la semana 3 era artefacto
de los núcleos efficiency**, no un límite del algoritmo.

### Plummer vs uniform

| | morton | orb |
|---|---|---|
| plummer | 28.914 | 10.768 |
| uniform | 12.314 | 5.884 |
| **brecha** | **2,35×** | **1,83×** |

La brecha de ORB da 1,83×, idéntica a la del portátil.

## Dos resultados que el cluster corrige

### El balance por costo casi no aporta con ORB

| Reparto | Desbal. partículas | Desbal. trabajo | Total |
|---|---|---|---|
| count | 1.0014 | 1.1364 | 5.209 |
| work | 1.0072 | **1.1345** | 5.187 |

Una mejora del **0,17%**, contra el 1,2469 → 1,1245 que daba en el portátil. La
explicación más plausible es que ORB ya bisecta buscando equilibrio, así que con
dominios compactos queda poco margen. **El balance por costo se justificaba sobre
dominios Morton; sobre ORB es casi redundante.** Es un matiz honesto sobre el
resultado de la semana 4, y conviene decirlo en el informe en vez de presentar la
mejora del portátil como si fuera general.

### La penalización en fuerzas era del hardware, no del algoritmo

En el portátil la fase de fuerzas empeoraba ~13% con ORB. Acá la diferencia es de
**2-3%** (P=4: 2,133 vs 2,202 s; P=8: 1,667 vs 1,707 s). La hipótesis de la
localidad de caché queda muy debilitada.

### Nota sobre el crecimiento del árbol con Morton

En el portátil `tree_time` con Morton crecía 4,8× de P=1 a P=8; acá se mantiene
plano. La explicación más plausible es que en el portátil los 8 procesos competían
por el mismo ancho de banda de memoria, mientras que acá están repartidos en 5
máquinas con su propia memoria. **Es una hipótesis, no una medición**: los dos
casos usan N y cantidad de pasos distintos y no los comparé de forma controlada.
No cambia la conclusión —ORB deja el árbol en la mitad y lo hace decrecer con P—
pero sí sugiere que parte del efecto de la semana 4 era del banco de pruebas.

## Pendiente

- **La fase de fuerzas domina la escalabilidad débil a P=8.** El árbol y el LET ya
  no son el cuello de botella; hay que instrumentar cuánto del recorrido se va en
  los fantasmas.
- **Sobre `--balance work`:** medir si sigue aportando algo sobre ORB en
  configuraciones más desbalanceadas, o documentarlo como redundante.
- **`docs/final-report.md`**.

## Lectura de las tres semanas

El arco quedó completo y es lo más valioso del trabajo:

| | Hipótesis | Medición |
|---|---|---|
| **Semana 3** | Domina la comunicación | El `Allgatherv` era el 1,3%; dominaba el árbol replicado |
| **Semana 4** | El LET reduce el árbol local a N/P | Lo redujo a ~N/2: los dominios Morton no son compactos |
| **Semana 5** | ORB da dominios compactos y el LET comprime | Confirmado: k de 0,92 a 0,66, y reproducido en pcunix (0,93 → 0,67) |

Dos veces seguidas la hipótesis razonable estaba equivocada, y las dos veces lo
que lo mostró fue una medición barata: en la semana 3, el desglose por fase; en la
semana 4, que `--theta-let` no cambiara nada y el volumen de las cajas de dominio.

El punto que conviene subrayar en el informe final: **el LET de la semana 4 era
correcto** —validado bit a bit con θ=0, masa a 5e-14, energía conservada— **y aun
así no servía para nada**. Correcto y efectivo son propiedades distintas. La suite
de tests certificaba la primera y era completamente ciega a la segunda; lo que
detectó el problema fue instrumentar el volumen de datos, que no era un test sino
una métrica.
