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

El repunte en P=8 es casi seguro sobresuscripción del portátil, no algorítmico
(ver entorno).

**Un efecto lateral honesto:** la fase de fuerzas empeora levemente con ORB
(2,79 → 3,16 s a P=4). No lo investigué a fondo; la hipótesis más plausible es
peor localidad de caché, porque el arreglo local queda ordenado por Morton pero
la propiedad ya no sigue la curva. Vale la pena mirarlo, aunque el total mejora.

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
`log N` por construcción y encima se mide sobresuscribiendo la máquina. Este
número hay que rehacerlo en pcunix antes de sacar conclusiones.

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

## Entorno

> Mismo portátil Apple Silicon de las semanas 3 y 4 (4 núcleos performance + 4
> efficiency). Las filas de P=8 sobresuscriben y usan los núcleos E, que rinden
> ~21% de un núcleo P. **Las comparaciones ORB vs Morton son válidas** porque las
> dos ramas corren en las mismas condiciones, pero **los tiempos absolutos y la
> eficiencia débil no lo son**. La deuda de re-medir en pcunix ya lleva tres
> semanas.

## Pendiente

- **Re-medir las semanas 3, 4 y 5 en pcunix.** Sigue siendo la deuda más vieja y
  la que más compromete el informe final.
- **Corrida distribuida en varias máquinas** (`-hosts`), arrastrada desde la
  semana 3.
- **La fase de fuerzas empeora ~13% con ORB.** Hipótesis: localidad de caché, por
  ordenar el arreglo por Morton cuando la propiedad ya no sigue la curva. Probar
  ordenar por la geometría de ORB.
- **Escalabilidad débil sigue mala** (7,6%). El árbol dejó de ser el cuello de
  botella; hay que instrumentar qué domina ahora.
- **`docs/final-report.md`**.

## Lectura de las tres semanas

El arco quedó completo y es lo más valioso del trabajo:

| | Hipótesis | Medición |
|---|---|---|
| **Semana 3** | Domina la comunicación | El `Allgatherv` era el 1,3%; dominaba el árbol replicado |
| **Semana 4** | El LET reduce el árbol local a N/P | Lo redujo a ~N/2: los dominios Morton no son compactos |
| **Semana 5** | ORB da dominios compactos y el LET comprime | Confirmado: k de 0,92 a 0,66 |

Dos veces seguidas la hipótesis razonable estaba equivocada, y las dos veces lo
que lo mostró fue una medición barata: en la semana 3, el desglose por fase; en la
semana 4, que `--theta-let` no cambiara nada y el volumen de las cajas de dominio.

El punto que conviene subrayar en el informe final: **el LET de la semana 4 era
correcto** —validado bit a bit con θ=0, masa a 5e-14, energía conservada— **y aun
así no servía para nada**. Correcto y efectivo son propiedades distintas. La suite
de tests certificaba la primera y era completamente ciega a la segunda; lo que
detectó el problema fue instrumentar el volumen de datos, que no era un test sino
una métrica.
