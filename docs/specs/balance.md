# Spec: Balance dinámico por costo

Etapa 9. Extiende `domain_partition` para repartir **trabajo** en vez de
**partículas**.

## Por qué

La semana 3 midió que los dos desbalances difieren:

| Métrica | Valor |
|---|---|
| Desbalance de **partículas** (max/avg) | 1.0000 |
| Desbalance de **trabajo** (max/avg) | 1.1149 |

`domain_partition` reparte los conteos casi perfectamente, pero el costo del
recorrido Barnes-Hut no es uniforme: una partícula en zona densa abre muchos más
nodos que una en la periferia. Repartir por conteo deja al proceso que recibe el
core haciendo ~11% más trabajo que el promedio.

## Medición del costo

`accumulate_force` lleva un contador de interacciones evaluadas (partículas de
hoja visitadas + nodos aceptados) y `compute_forces_bh_range` lo guarda en
`p[i].work`. Un incremento por interacción es despreciable frente al `sqrt` que ya
hay en cada una.

El peso vive **dentro** de `Particle`, no en un arreglo paralelo: así viaja solo
en la migración, sin que haya que mantener sincronizada una segunda permutación.
Cuesta 8 bytes por partícula.

> ⚠ Agregar un campo al struct **obliga a actualizar `mpi_particle_type_init`**.
> Si no se hace, la migración copia basura en silencio. El test 11 (checksum de
> identidad) lo detecta.

## Algoritmo

`domain_partition_ex(d, p, n_local, n_global, use_work, comm)`.

No cambia nada del esquema: es la misma bisección simultánea de los P-1
splitters, con la misma cantidad de colectivos (~63). Solo cambia la magnitud que
se acumula — de "cuántas partículas hay debajo de esta clave" a "cuánto trabajo
hay debajo de esta clave".

Como el arreglo local ya está ordenado por clave, alcanza con un **prefijo de
`work`** calculado una vez por llamada:

```
prefix[0] = 0
prefix[i+1] = prefix[i] + p[i].work
trabajo con clave < K  =  prefix[count_below(K)]
```

Mismo O(log n) por candidato que el conteo, con una lectura extra.

## Casos borde

- **Primer paso:** `work` vale 0 para todas. Se detecta `w_global <= 0` y se cae
  al reparto por conteo, en vez de dividir por cero o dejar los splitters en 0.
- **Peso retrasado un paso:** el costo es el del paso anterior. Es la práctica
  estándar y es válida porque la distribución cambia poco entre pasos
  consecutivos. Conviene decirlo explícitamente en el informe.

## Resultado esperado (y contraintuitivo)

**El desbalance de partículas debe empeorar.** Los procesos que reciben zonas
densas se quedan con menos partículas, que es precisamente el punto. La métrica
de éxito es el desbalance de trabajo, no el de conteo.

Medido con N=100K, LET:

| P | Reparto | Desbal. partículas | Desbal. trabajo |
|---|---|---|---|
| 4 | count | 1.0015 | 1.1382 |
| 4 | **work** | 1.0723 | **1.0473** |
| 8 | count | 1.0030 | 1.2469 |
| 8 | **work** | 1.2142 | **1.1245** |

Con reparticiones más frecuentes (test 17): **1.2072 → 1.0013**.

A P=8 la mejora se traduce en 7% menos tiempo total.

## CLI

```
--balance work|count       # default: work
--rebalance-every K        # cada cuántos pasos se recalculan los splitters
```

Los dos modos conviven a propósito: sin poder correr ambos no hay forma de medir
la mejora, y la comparación es uno de los resultados del informe.

## Reutilizable para ORB

La bisección por histograma sirve igual para bisectar una **coordenada** en vez
del espacio de claves, que es lo que haría falta para ORB — el trabajo futuro
prioritario según [let.md](let.md).
