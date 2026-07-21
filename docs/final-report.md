# Simulación paralela del problema de N cuerpos con MPI + OpenMP

## Reporte final — base para el informe del proyecto

---

## 1. El problema y por qué necesita paralelismo

Simular la evolución gravitatoria de N cuerpos en 3D. Cada partícula siente la
atracción de todas las demás, de modo que la evaluación directa de fuerzas cuesta
**O(N²)** por paso: con N = 100.000 son 10¹⁰ interacciones por paso, y una
simulación útil necesita miles de pasos.

Hay dos ejes independientes para atacarlo, y el trabajo usa los dos:

| Eje | Herramienta | Qué reduce |
|---|---|---|
| **Algorítmico** | Barnes-Hut | El costo por paso, de O(N²) a O(N log N) |
| **Paralelo** | MPI + OpenMP | El tiempo de pared, repartiendo ese costo |

Los dos son necesarios y **el algorítmico va primero**: medido en este trabajo,
Barnes-Hut da **14,8× sobre el método directo con N = 100.000** en un solo núcleo.
Ningún reparto entre 8 procesos habría conseguido eso. Paralelizar un algoritmo
malo es la forma más cara de no resolver un problema.

---

## 2. La solución construida

### 2.1 El núcleo secuencial

**Integrador leapfrog en forma KDK** (kick–drift–kick). Es simpléctico: conserva
la energía a largo plazo en vez de derivar. Medido: drift relativo de energía de
**9,0e-06** en corridas de 50.000 partículas.

**Octree y criterio de apertura.** El espacio se subdivide recursivamente en
octantes hasta que cada hoja tiene a lo sumo una partícula. Cada nodo guarda la
masa total y el centro de masa de su subárbol. Al calcular la fuerza sobre una
partícula se recorre el árbol desde la raíz: si una celda de lado `s` está a
distancia `d` y cumple **s/d < θ**, todo su contenido se aproxima por un único
cuerpo en su centro de masa; si no, se desciende a los 8 hijos.

El árbol usa un **pool de nodos con índices enteros** en lugar de punteros. Reduce
las llamadas a `malloc`, mejora la localidad de caché y —lo que importa acá— es
directamente serializable para enviarlo por MPI.

**Propiedad clave:** con **θ = 0** no se aproxima nada y el resultado coincide con
el método directo. Esta propiedad, que parece un caso degenerado sin interés,
terminó siendo la herramienta de validación más importante del proyecto (§4).

### 2.2 El reparto híbrido: qué hace cada capa

```
        ┌─────────────── MPI: reparte PARTÍCULAS entre procesos ───────────────┐
        │                                                                      │
   proceso 0          proceso 1          proceso 2          proceso 3
   caja del           caja del           caja del           caja del
   dominio            dominio            dominio            dominio
        │                  │                  │                  │
   ┌────┴────┐        ┌────┴────┐        ┌────┴────┐        ┌────┴────┐
   h0   h1  ...       h0   h1  ...       ...                ...
   └── OpenMP: reparte el RECORRIDO del árbol entre hilos ──┘
```

**Por qué esta división y no otra.** Los dos niveles resuelven cosas distintas:

- **MPI** ataca el problema de que los datos no entran, o no conviene que entren,
  en una sola máquina. Cada proceso es dueño de una región del espacio y solo
  guarda sus partículas. Es lo que permite crecer más allá de un nodo.
- **OpenMP** ataca el hecho de que dentro de un proceso el recorrido del árbol es
  **embarazosamente paralelo**: cada partícula lee un árbol inmutable y escribe
  solo su propia aceleración. No hay carreras ni acumulación compartida, y los
  hilos comparten el árbol sin duplicarlo.

Esa asimetría —los hilos comparten el árbol, los procesos tienen que
intercambiarlo— es el eje de casi todos los hallazgos del trabajo.

### 2.3 Descomposición del dominio: ORB

Cada proceso es dueño de una **caja** del espacio, obtenida por **bisección
recursiva ortogonal**: en cada nivel se corta por la dimensión más larga de la
caja, buscando el plano que reparte el trabajo por la mitad entre las dos mitades
del grupo de procesos.

```
orb_build(caja, procesos [first, first+n)):
    si n == 1: hoja, first es el dueño de la caja
    eje = dimensión más larga
    bisecar la coordenada hasta repartir el trabajo nL/n
        (cada iteración: un MPI_Allreduce del trabajo por debajo del candidato)
    recurrir en las dos mitades
```

El árbol ORB tiene **2P−1 nodos y está replicado en todos los procesos**. Son
cientos de bytes, así que cada proceso conoce la caja de todos los demás **sin
comunicación**, que es exactamente lo que necesita el intercambio de árboles.

Dos detalles de implementación que evitaron bugs:

- **La pertenencia se lleva por el camino recorrido, no por un test geométrico.**
  El arreglo de índices locales se particiona in-place en cada nivel con el mismo
  predicado `<` que usa la función de dueño. Si el conteo de la bisección y la
  asignación usaran comparaciones distintas, una partícula justo sobre el plano se
  contaría de un lado y se enviaría al otro, y se perdería. Llevarlo por el camino
  lo vuelve **imposible por construcción** en vez de dependiente de que dos
  comparaciones coincidan.
- **La función de dueño no hace test de contención**: cualquier posición, aun
  fuera de la caja global, obtiene un dueño válido. Es lo que hace seguro el caso
  de una partícula que se escapa del dominio entre reparticiones.

> **Nota histórica importante.** Las primeras versiones repartían por **rangos
> contiguos de clave Morton** (curva Z). Funciona para migrar y es más simple,
> pero produce dominios que no son compactos, y eso hundió el rendimiento por
> razones que tardaron una semana entera en aparecer. Es el Hallazgo 3 (§5.3).

### 2.4 Migración

Después de mover las partículas, algunas quedan fuera de la caja de su proceso.
La migración las reubica:

```
1. dueño de cada partícula = descenso por el árbol ORB según su posición
2. agrupar el arreglo local por destino (counting sort estable)
3. MPI_Alltoall     -> cuántas voy a recibir de cada uno
4. MPI_Alltoallv    -> el intercambio
5. reordenar por clave Morton (ahora solo por localidad de caché)
```

### 2.5 LET: el árbol localmente esencial

El problema central del N-cuerpos distribuido: **para calcular sus fuerzas, cada
proceso necesita información de partículas que no tiene.**

La solución ingenua —que este trabajo implementó primero, a propósito, como
andamio— es que cada proceso replique el arreglo global con un `MPI_Allgatherv` y
construya el árbol completo. Funciona, es trivial de validar, y **no escala**.

El LET manda, en cambio, solo el subconjunto del árbol que cada destino necesita
para calcular sus fuerzas **con la misma precisión** que tendría con el árbol
entero. Lo lejano llega resumido en unos pocos nodos agregados.

**Criterio de exportación (MAC conservador).** Para el nodo `c` de lado `s` y
centro de masa `cm`, y la caja destino `B_r`:

```
d_min = distancia mínima de cm a B_r      (0 si cm cae dentro de B_r)
exportar el resumen si   s² < θ² · (d_min² + ε²)
```

**Por qué es correcto:** el receptor, al calcular la fuerza sobre su partícula
`i ∈ B_r`, aplicaría `s² < θ²·(d_i² + ε²)` con `d_i ≥ d_min`. Usar la distancia
**mínima** sobre toda la caja garantiza que nunca se exporta un resumen donde el
receptor quería detalle. El caso inverso —abrir un nodo que el receptor habría
aceptado— manda datos de más: cuesta ancho de banda, no precisión.

**Representación.** Lo exportado viaja como una `Particle` común: los nodos
aceptados como pseudo-partículas (`pos = cm`, `mass` = masa del subárbol,
`id = −1` como marca de fantasma), las hojas como partículas reales. Así se
reutilizan el tipo MPI derivado y `MPI_Alltoallv` **sin ningún empaquetado a
mano**.

**Los fantasmas son transitorios por diseño.** Se agregan al final del arreglo
local y se insertan en el árbol ya construido —no se reconstruye nada, porque la
caja raíz es la misma—; las fuerzas se calculan solo para el rango local, y la
migración del paso siguiente los descarta sola. Todo lo que sea checksum, energía
o migración usa el rango de las partículas propias.

### 2.6 Balance de carga por costo

El desbalance de **partículas** y el de **trabajo** no son lo mismo: una partícula
en una zona densa abre muchos más nodos del árbol que una en la periferia. Medido
en la versión replicada: desbalance de partículas **1,0000** contra desbalance de
trabajo **1,1149**.

La solución: el recorrido cuenta las interacciones que evalúa y las guarda en la
propia partícula; la bisección de ORB reparte **esa** magnitud en vez del conteo.
El peso vive dentro del struct y no en un arreglo paralelo, para que viaje solo en
la migración sin mantener una segunda permutación sincronizada.

**Señal esperada y contraintuitiva:** al activar el reparto por costo, el
desbalance de *partículas* empeora mientras el de *trabajo* mejora. Es el objetivo,
no un problema.

### 2.7 El paso completo

```
kick(dt/2)  ·  drift(dt)                             [OpenMP]

[cada K pasos] bbox global (Allreduce MIN/MAX)  ->  ORB (Allreduce × bisección)
migración                                             [Alltoall + Alltoallv]

árbol local sobre la caja raíz compartida
centros de masa
intercambio LET                                       [Allgather + Alltoallv]
fantasmas insertados en el mismo árbol
centros de masa (cambiaron)

fuerzas sobre el rango propio                         [OpenMP]

kick(dt/2)                                            [OpenMP]
```

**El orden importa:** migrar *antes* de construir el árbol. Si se migrara al final
del paso, el árbol se construiría con partículas que ya no pertenecen a nadie.

---

## 3. Uso de MPI y OpenMP: decisiones y justificación

### 3.1 MPI — solo colectivos, ningún punto a punto

| Primitiva | Uso |
|---|---|
| `MPI_Type_create_struct` + `_resized` | Tipo derivado para `Particle` |
| `MPI_Allreduce` (MIN/MAX/SUM) | Bbox global, bisección ORB, energías, checksums |
| `MPI_Allgather` | Cajas de dominio para el LET |
| `MPI_Alltoall` | Intercambio de conteos previo a los `v` |
| `MPI_Alltoallv` | Migración y LET |
| `MPI_Scatterv` | Reparto inicial |
| `MPI_Allgatherv` | Réplica global (esquema de referencia) |
| `MPI_Barrier` / `MPI_Wtime` | Instrumentación por fase |

**Decisión: todo colectivo, cero `Send`/`Recv`.** Un patrón punto a punto armado a
mano para la migración habría sido más código y, sobre todo, **puede deadlockear
por un orden de envíos mal elegido**. `MPI_Alltoallv` no puede: es colectivo, y
además el runtime elige internamente el algoritmo (bruck para mensajes chicos,
pairwise para grandes) según el tamaño. **En todo el proyecto no hubo un solo
deadlock.**

**Decisión: tipo derivado, nunca `MPI_BYTE` sobre el struct crudo.** Enviar el
struct como bytes funciona solo si todos los nodos son homogéneos, y esconde
errores de alineación. Se describen los campos con `MPI_Type_create_struct`
midiendo los desplazamientos sobre una instancia real (para que el padding que
elija el compilador quede descrito), y se fija el extent con
`MPI_Type_create_resized` a `sizeof(Particle)`, que es imprescindible para enviar
arreglos.

> **Trampa concreta encontrada:** agregar un campo al struct obliga a actualizar el
> tipo derivado. Si no se hace, la migración **copia basura en silencio**. Es el
> error más fácil de cometer y el más difícil de ver. Lo detecta el test de
> checksum de identidad.

**Decisión: bisección simultánea.** El particionamiento busca P−1 cortes a la vez
en un único `Allreduce` de P−1 enteros por iteración, en lugar de P−1 reducciones
independientes. El costo pasa de ~63·(P−1) colectivos a ~63.

### 3.2 OpenMP — dónde sí y dónde no

| Región | Directiva | Razón |
|---|---|---|
| Recorrido Barnes-Hut | `parallel for schedule(dynamic, 64)` | Costo por partícula muy desigual |
| Energía potencial | `... reduction(+:U)` | Acumulador compartido |
| `kick` / `drift` | `parallel for schedule(static)` | Costo uniforme |
| Construcción del octree | **no paralelizada** | Inserción con `realloc` sobre un pool compartido |

**Por qué `dynamic` y no `static` en el recorrido:** el costo por partícula no es
uniforme —una partícula en zona densa abre muchos más nodos—, así que con reparto
estático los hilos que tocan zonas vacías terminan primero y quedan ociosos. El
chunk de 64 amortiza el costo de sincronización del planificador.

**Determinismo como propiedad de diseño.** Cada iteración escribe **solo** su
propia aceleración, y el orden de las sumas dentro de una partícula no depende del
número de hilos. Por lo tanto el resultado es **bit a bit idéntico** al secuencial
para cualquier cantidad de hilos. Eso no es una curiosidad: convierte "¿hay una
carrera de datos?" en un test automatizable con `memcmp` (§4.2). Las reducciones
de energía sí dependen del orden, así que ahí se usa tolerancia y no igualdad
exacta.

**Lo que no se paralelizó, y por qué se midió aparte:** la construcción del octree
es inserción secuencial en un pool que se reasigna con `realloc`. Paralelizarla
exigiría una construcción bottom-up por claves Morton, fuera del alcance. Se dejó
serial **y se instrumentó como fase propia**, porque es el techo de Amdahl de la
capa OpenMP.

---

## 4. Metodología de validación

Es, junto con los hallazgos, la parte más transferible del trabajo. **19 tests**,
construidos en capas, de propiedades físicas a comparaciones contra oráculos.

### 4.1 La escalera de oráculos

Cada nivel valida al siguiente:

```
propiedades analíticas (órbita de 2 cuerpos, conservación de momento)
        ↓
método directo O(N²)         ← oráculo de Barnes-Hut
        ↓
Barnes-Hut secuencial        ← oráculo de la versión distribuida
        ↓
réplica global distribuida   ← oráculo del LET
        ↓
LET
```

### 4.2 Dos tests de igualdad bit a bit, y por qué valen tanto

**Determinismo de OpenMP.** Calcular las aceleraciones con 1, 2, 4 y 8 hilos y
exigir `memcmp` idéntico. Cualquier diferencia, por chica que sea, delata una
carrera o una acumulación compartida. Es infinitamente más sensible que comparar
con tolerancia.

**LET con θ = 0.** Con θ = 0 el criterio de exportación no acepta ningún resumen,
el LET degenera en réplica completa, el árbol resultante tiene la misma estructura
y el mismo orden de sumas que el secuencial, y la igualdad vuelve a ser **exacta**.

Este segundo test resultó ser **el guardián de todo el proyecto**. Detecta, con un
único número que debe ser exactamente cero:

- cualquier interacción que la selección del LET se saltee,
- cualquier partícula perdida o duplicada en la migración,
- cualquier desalineación de la caja raíz entre procesos.

Fue lo que dio confianza para reescribir por completo la descomposición del
dominio y la migración en la última etapa: mientras ese test siguiera dando cero,
la reescritura no había roto nada.

### 4.3 Invariantes globales baratos

- **Checksums Σid y Σid².** El segundo momento es necesario: Σid sola no detecta
  que dos partículas intercambien identidad.
- **Conservación de masa del LET.** La suma de masas que "ve" cada proceso
  (propias + fantasmas) debe dar la masa total del sistema, porque el LET resume
  pero no descarta. Medido: **4,9e-14**. Detecta casi cualquier bug de selección
  sin necesitar una referencia secuencial.
- **Conservación de partículas.** Σ n_local == N en todo paso.

### 4.4 Dos bugs que encontró la validación y no la inspección

**Codificación Morton sin acotar.** La coordenada normalizada se convertía a
entero sin acotarla. Con el bounding box congelado entre reparticiones, una
partícula que se sale del box producía: si excedía el máximo, el enmascarado
**envolvía** el valor y una partícula del extremo derecho recibía una clave del
extremo izquierdo, migrando al proceso equivocado; si quedaba por debajo del
mínimo, el double negativo convertido a entero sin signo es **comportamiento
indefinido**. No se manifestaba antes porque el box se recalculaba cada paso.

**Carrera de datos en la fuerza directa.** La versión O(N²) aprovecha la tercera
ley de Newton acumulando en `p[j].acc` dentro del bucle interno. Con
`parallel for` sobre `i`, dos hilos escriben el mismo `p[j].acc`: **carrera
silenciosa**. La solución fue una variante sin simetría (cada `i` escribe solo lo
suyo), al doble de operaciones pero reproducible bit a bit, conservando la versión
simétrica intacta como oráculo secuencial.

Los dos son del tipo que no se ve leyendo el código y que en producción se
manifiesta como "resultados raros de vez en cuando".

---

## 5. Hallazgos

Esta es la parte central del trabajo. **Tres veces la hipótesis razonable estaba
equivocada, y las tres veces lo mostró una medición barata.**

### 5.1 Hallazgo 1 — El cuello de botella no era la comunicación

La primera versión distribuida replicaba el arreglo global con `Allgatherv`. La
hipótesis de trabajo, y la del plan, era que ese `Allgatherv` dominaría el tiempo
y que ese sería el argumento para implementar el LET.

**La medición lo refutó.** Con 8 procesos:

| Fase | Tiempo | % del total |
|---|---|---|
| `Allgatherv` (comunicación) | 0,071 s | **1,3%** |
| Construcción del árbol | 1,133 s | 21% |
| Fuerzas | 3,270 s | 60% |

El costo real de la replicación no era mover los datos: era que **cada proceso
construía el árbol completo sobre las N partículas**. Trabajo idéntico, hecho P
veces en paralelo, que no se reparte entre nadie. `tree_time` crecía **4,8× de
P=1 a P=8**, y por eso el rendimiento **empeoraba** al pasar de 4 a 8 procesos.

**La lección, que es general:** en un código paralelo, el trabajo redundante puede
costar mucho más que la comunicación. Y no hay forma de saber cuál domina sin
**instrumentar por fase**. La intuición de que "en distribuido lo caro es la red"
es un prejuicio que costó una semana.

### 5.2 Hallazgo 2 — Correcto y efectivo son propiedades distintas

El LET se implementó y **validó completamente**:

| Verificación | Resultado |
|---|---|
| θ = 0 contra el secuencial | **bit a bit exacto** |
| Conservación de masa | 4,9e-14 |
| Drift de energía, 200 pasos | 1,21e-04 |
| Suite de tests | 17/18 |

Y sin embargo **no servía para nada**. El volumen de datos importados crecía como
**N^0,92** —prácticamente lineal, la misma asintótica que la réplica que venía a
reemplazar—. El objetivo entero de la etapa no se cumplía.

**La suite de tests era completamente ciega a esto.** Todos los tests medían
*corrección*, y el LET era correcto. Lo que detectó el problema fue una **métrica**
—el volumen de datos intercambiados— que no era un test y no tenía criterio de
aprobación.

**La lección:** en trabajo de rendimiento, hacen falta dos familias de
verificación separadas. Los tests dicen si el resultado está bien; **las métricas
dicen si el mecanismo está haciendo lo que se supone que hace.** Un test suite
verde puede convivir perfectamente con un componente que no cumple ninguna de sus
funciones.

### 5.3 Hallazgo 3 — La geometría de la descomposición determina el volumen de comunicación

¿Por qué el LET no comprimía? El diagnóstico salió de dos mediciones baratas:

**Primera pista:** subir el θ del criterio de exportación de 0,5 a 4,0 —un factor
8— **no cambiaba nada** en el volumen exportado. Si el criterio gobernara la
selección, eso tendría que haberlo desplomado. Que no se moviera significaba que
la rama que acepta el resumen prácticamente nunca se tomaba.

**Segunda:** instrumentar las cajas de dominio.

```
caja[1] lados=(  3,2   3,2   4,8)   volumen relativo = 0,0000
caja[3] lados=(402,5 293,2 261,3)   volumen relativo = 0,4456
```

La caja de un proceso cubría el **44,6% del volumen del dominio**. Contra un
destino así, la distancia mínima de casi cualquier nodo es cero, el criterio nunca
se cumple, y se exporta todo sin resumir.

**La causa raíz.** El supuesto equivocado era que **un tramo contiguo de clave
Morton describe una región compacta del espacio**. No es así. En una distribución
centralmente concentrada (Plummer), el pico de densidad cae justo donde se tocan
los octantes de primer nivel, así que las claves del núcleo denso se reparten por
**todo** el espacio de claves. Cualquier corte por percentil le da a cada proceso
una mezcla de núcleo y halo, y los dominios se interpenetran geométricamente.

**Verificación independiente del diagnóstico:** con una distribución uniforme, sin
pico central, el mismo caso bajaba de 61.056 a 22.029 fantasmas. La patología era
de la distribución combinada con dominios no compactos, no del código.

**La corrección** fue reemplazar los rangos de curva por ORB, que produce dominios
que **son cajas por construcción**. Y el resultado más elocuente:

> **El módulo del LET no necesitó un solo cambio de lógica.** El criterio de
> exportación siempre había sido correcto; lo que fallaba era la geometría que
> recibía.

| | Morton | ORB |
|---|---|---|
| Exponente de crecimiento del volumen | **N^0,93** | **N^0,67** |
| Brecha Plummer / uniforme | 2,35× | 1,83× |

El 0,67 está muy cerca del **2/3 teórico de un efecto de superficie**, que es lo
que un LET debe costar: cada proceso necesita detalle solo cerca de la frontera de
su dominio.

**La lección, que es la más transferible del proyecto:** en un método distribuido
basado en localidad espacial, **la calidad geométrica de la descomposición
determina el volumen de comunicación**, y puede hacerlo por órdenes de magnitud.
Una descomposición que reparte bien la *carga* puede ser pésima para la
*comunicación*. Las curvas de llenado de espacio (Morton, Hilbert) preservan
localidad **en promedio**, y eso no alcanza cuando el criterio depende de la
distancia mínima a la envolvente del dominio.

### 5.4 Hallazgo 4 — El entorno de medición puede inventar fenómenos

Los primeros cuatro conjuntos de mediciones se hicieron en un portátil Apple
Silicon con **4 núcleos de rendimiento + 4 de eficiencia**. Al re-medir en el
cluster de fing (5 nodos homogéneos, ejecución realmente distribuida), tres
conclusiones cambiaron:

| Conclusión en el portátil | En hardware homogéneo |
|---|---|
| Speedup OpenMP 3,94× con 8 hilos: "no cumple el criterio de 5×" | **1,96× con 2 hilos = 98% de eficiencia.** Era el hardware, no el algoritmo |
| ORB penaliza la fase de fuerzas un 13% (¿localidad de caché?) | **2-3%.** Era del banco de pruebas |
| Speedup MPI nunca supera 2,80× con P=4 | **3,20× con P=4** |
| Eficiencia débil 7,6% con P=8 | **62% con P=4**, 33% con P=8 |

**Lo que sí se sostuvo intacto:** el resultado central. El exponente de crecimiento
del volumen del LET dio **0,67 en el cluster contra 0,66 en el portátil**, y la
serie de cocientes coincidió hasta la segunda decimal en los cuatro tamaños, con
otro compilador (GCC 15.2 en vez de clang), otra implementación de MPI (MPICH en
vez de Open MPI) y otra arquitectura (x86-64 en vez de ARM).

**La lección:** las comparaciones *relativas* medidas en las mismas condiciones son
robustas; los *números absolutos* y todo lo que dependa del comportamiento de los
núcleos, no. Un banco de pruebas con núcleos heterogéneos o sobresuscripto puede
generar fenómenos que no existen, y hacer perder tiempo buscándoles explicación
algorítmica. **Conviene decir siempre en qué hardware se midió.**

### 5.5 Hallazgo 5 — El reparto por costo depende de la descomposición

El balance por costo estimado —repartir el trabajo medido en vez del número de
partículas— fue una mejora clara sobre dominios Morton:

| | Desbalance de partículas | Desbalance de trabajo |
|---|---|---|
| Reparto por conteo | 1,0014 | 1,2072 |
| Reparto por costo | 1,0640 | **1,0013** |

Pero sobre dominios ORB, medido en el cluster, la mejora fue del **0,17%**
(1,1364 → 1,1345). La explicación más plausible es que ORB ya bisecta buscando
equilibrio, así que sobre dominios compactos queda poco margen que ganar.

**La lección:** una optimización puede volverse redundante cuando cambia otra
parte del sistema. Conviene re-medir las mejoras viejas después de cada cambio
estructural, en lugar de asumir que se acumulan.

### 5.6 Hallazgo 6 — El punto óptimo del reparto proceso/hilo

Con la misma cantidad de núcleos físicos (8):

| Configuración | Tiempo | Desbalance de trabajo |
|---|---|---|
| 8 procesos × 1 hilo | 2,767 s | 1,1303 |
| **4 procesos × 2 hilos** | **1,871 s** | **1,0047** |

**Menos procesos con más hilos gana por 1,48×**, y la causa está en la última
columna: con dominios más grandes el reparto es más parejo. Además, cada proceso
adicional agrega un dominio más al que exportar, y por lo tanto más comunicación.

La regla práctica que se desprende: **usar OpenMP para llenar los núcleos dentro
de un nodo y MPI para cruzar nodos**, en vez de lanzar un proceso MPI por núcleo.

---

## 6. Resultados consolidados

Medidos en el cluster de fing: 5 nodos Intel i3-4170 (2 núcleos físicos c/u), GCC
15.2.1, MPICH, x86-64. Los barridos van distribuidos a 2 procesos por nodo, de
modo que **P = 1..8 son 1 a 8 núcleos físicos sin sobresuscribir**.

### Escalabilidad fuerte (N = 50.000)

| P | Morton | speedup | ORB | speedup | eficiencia |
|---|---|---|---|---|---|
| 1 | 10,142 s | 1,00× | 10,118 s | 1,00× | 100% |
| 2 | 5,325 | 1,90× | 5,190 | 1,95× | 97% |
| 4 | 4,371 | 2,32× | 3,161 | **3,20×** | 80% |
| 8 | 4,176 | 2,43× | 2,777 | 3,64× | 46% |

### Escalabilidad débil (N/P = 12.500)

| P | Morton | ORB | eficiencia ORB |
|---|---|---|---|
| 2 | 2,343 s | 2,294 s | 85% |
| 4 | 4,526 | 3,164 | **62%** |
| 8 | 8,885 | 5,980 | 33% |

### OpenMP (2 núcleos físicos)

| Hilos | Speedup total | Speedup de la fase de fuerzas |
|---|---|---|
| 2 | **1,96×** | 2,00× |
| 4 (hyperthreading) | 2,39× | 2,43× |

### Volumen del LET (P = 4)

| N | Morton | cociente | ORB | cociente |
|---|---|---|---|---|
| 12.500 | 8.544 | 2,73 | 4.225 | **1,35** |
| 25.000 | 17.102 | 2,74 | 6.804 | **1,09** |
| 50.000 | 28.914 | 2,31 | 10.768 | **0,86** |
| 100.000 | 58.739 | 2,35 | 16.906 | **0,68** |

A N fijo subiendo P, los fantasmas de ORB se mantienen **casi constantes** (9.038 /
10.769 / 9.794 para P = 2/4/8): la firma de un término de superficie.

> Un número elocuente: con Morton y P=2, cada proceso importa **24.995** fantasmas
> sobre N=50.000. Es decir, **exactamente todas las partículas del otro proceso**.
> Compresión cero.

### Corrección

**19/19 tests** con P = 2, 3, 4 y 8 distribuidos en 5 máquinas, más 9/9 en la
versión secuencial. Sin fugas de memoria. Sin deadlocks. Idénticos resultados en
dos arquitecturas, dos compiladores y dos implementaciones de MPI, **incluidos los
dos tests de igualdad bit a bit**.

---

## 7. Lecciones metodológicas para trabajo con MPI + OpenMP

1. **Instrumentar por fase antes de optimizar.** Sin el desglose
   árbol/fuerzas/comunicación/migración, el diagnóstico habría sido el equivocado
   dos veces. La intuición sobre dónde está el costo no es confiable.

2. **Separar tests de métricas.** Los tests verifican corrección; las métricas
   verifican que el mecanismo cumpla su función. Una suite verde es perfectamente
   compatible con un componente inútil.

3. **Buscar propiedades que permitan igualdad bit a bit.** El determinismo de
   OpenMP y el caso θ = 0 del LET convirtieron dos preguntas difíciles —"¿hay una
   carrera?", "¿falta alguna interacción?"— en tests automáticos con un criterio
   sin ambigüedad. Vale la pena diseñar buscando esas propiedades.

4. **Conservar el camino anterior detrás de un flag.** Cada esquema reemplazado
   quedó disponible por línea de comandos (`--decomp`, `--exchange`, `--balance`).
   **No es deuda técnica: es el instrumento de medición.** Sin poder correr las dos
   versiones bajo condiciones idénticas, ninguna de las comparaciones de este
   informe existiría.

5. **Preferir colectivos.** Cero deadlocks en todo el proyecto, y el runtime elige
   mejores algoritmos internos que un patrón armado a mano.

6. **Describir los tipos, no enviar bytes.** Y recordar que agregar un campo a un
   struct obliga a tocar el tipo derivado, o se copia basura en silencio.

7. **Documentar el hardware de medición.** Un banco de pruebas heterogéneo o
   sobresuscripto inventa fenómenos.

8. **En métodos basados en localidad espacial, la geometría de la descomposición
   es una decisión de primer orden**, no un detalle de implementación.

---

## 8. Limitaciones y trabajo futuro

- **La fase de fuerzas domina la escalabilidad débil con P = 8** (33%). El árbol y
  el LET dejaron de ser el cuello de botella; falta instrumentar cuánto del
  recorrido se consume en los nodos fantasma importados.
- **La construcción del octree no está paralelizada** con OpenMP. Es el techo de
  Amdahl de la capa de hilos. Requeriría construcción bottom-up por claves Morton.
- **El reparto por costo es casi redundante sobre ORB.** Falta medir si aporta algo
  en configuraciones más desbalanceadas, o documentarlo como innecesario.
- **Un solo tipo de distribución inicial evaluado en profundidad** (Plummer, más
  comprobaciones puntuales con cubo uniforme). Un colapso frío o una distribución
  con subestructura ejercitarían más el rebalanceo dinámico.
- **Escala máxima probada: 8 núcleos físicos.** Las conclusiones sobre el
  comportamiento asintótico se apoyan en el exponente de crecimiento medido, no en
  corridas a gran escala.

---

## 9. Conclusión

Se construyó un simulador de N cuerpos híbrido MPI + OpenMP con Barnes-Hut,
descomposición del dominio por bisección recursiva ortogonal, intercambio de
árboles localmente esenciales y balance dinámico por costo, validado con 19 tests
en dos arquitecturas.

Pero el resultado más valioso no es el código, sino **el patrón que se repitió tres
veces**: una hipótesis razonable sobre dónde estaba el costo, sostenida por
intuición y no por datos, refutada por una medición barata.

- Se esperaba que dominara la comunicación → dominaba el cómputo redundante.
- Se esperaba que el LET redujera el árbol local a N/P → lo redujo a ~N/2, porque
  los dominios no eran compactos.
- Se esperaba que el hardware homogéneo confirmara los límites medidos → mostró
  que varios eran artefactos del banco de pruebas.

En ninguno de los tres casos el problema se descubrió leyendo código o razonando
sobre el algoritmo. En los tres, lo mostró instrumentar algo y mirar el número.
**Esa es la conclusión práctica del trabajo sobre programación paralela: el ciclo
medir-diagnosticar-corregir no es una fase final de ajuste, es el método.**
