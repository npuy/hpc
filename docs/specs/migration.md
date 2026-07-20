# Módulo: migration

## Descripción general

Reubica entre procesos las partículas que, al moverse, salieron del tramo de claves Morton de su dueño (etapa 7 del plan). Contiene además `replicate_particles`, la réplica global que la semana 3 usa **en lugar del LET**, y `particle_checksum`, la verificación de que ninguna partícula se pierde.

## Conceptos

### Por qué `Alltoallv` y no `Send`/`Recv`

El patrón de comunicación es irregular: cada proceso envía una cantidad distinta de partículas a cada otro, y esa cantidad cambia en cada paso. Armarlo a mano con `MPI_Send`/`MPI_Recv` exige elegir un orden de envíos que no se trabe, y un orden mal elegido produce un deadlock que aparece de forma intermitente y solo con ciertos `P`.

`MPI_Alltoallv` es colectivo: no puede deadlockear por orden de envíos, y el runtime elige internamente el algoritmo (bruck para mensajes chicos, pairwise para grandes) según el tamaño. Es menos código y menos riesgo, sin costo de rendimiento a esta escala.

### Orden de las fases dentro del paso

La migración va **antes** de construir el árbol, no al final del paso:

```
kick(dt/2) → drift(dt) → [reparticionar cada K] → migrar → replicar → arbol → fuerzas → kick(dt/2)
```

Si se migrara al final, el árbol de ese paso se construiría con partículas que ya no pertenecen a su proceso, y el rango local `[offset, offset+n_local)` dentro del arreglo replicado dejaría de ser contiguo — que es justo la propiedad de la que depende `compute_forces_bh_range`.

### Algoritmo de `migrate_particles`

```
1. recodificar claves Morton con el bounding box global congelado
2. dest[i] = domain_owner(clave[i])           y contar send_counts[P]
3. desplazamientos por prefijo               -> send_displs[P]
4. agrupar por destino (counting sort estable sobre dest[])
5. MPI_Alltoall(send_counts -> recv_counts)   -- cuanto recibe cada uno
6. MPI_Alltoallv(sendbuf -> recvbuf)
7. reemplazar el arreglo local y reordenar por clave
```

El paso 7 restaura la invariante de orden que `domain_partition` asume. Recodificar las claves con el mismo box congelado es idempotente (el intercambio no cambia posiciones), así que ese `morton_sort` solo reordena.

### La invariante que hace innecesario un sort global

Cada proceso mantiene sus partículas ordenadas por clave, y los tramos de claves son disjuntos y crecientes con el rank. Por lo tanto **la concatenación que produce `Allgatherv` ya está globalmente ordenada por Morton**: no hace falta ningún sort posterior, y las partículas locales ocupan el rango contiguo `[offset, offset + n_local)`.

### `replicate_particles` es provisional

Cada proceso replica el arreglo global completo (`MPI_Allgatherv`) y construye el árbol entero, pero solo calcula fuerzas para su tramo. Esto hace que las fuerzas sean **numéricamente idénticas** a la versión secuencial — el test 13 mide error exactamente `0.0`.

El costo es O(N) en memoria y comunicación por proceso, y sobre todo la **construcción del árbol duplicada P veces**. Medido en el proyecto, esa duplicación —y no el `Allgatherv`, que resulta ser barato— es lo que limita la escalabilidad a ~4 procesos (ver `docs/week3-report.md`).

Está aislada en una sola función a propósito: la semana 4 la sustituye por el intercambio de árboles localmente esenciales tocando únicamente este punto.

### Checksums de identidad

`particle_checksum` calcula `Σ id` y `Σ id²` globales. Si ambos coinciden con sus valores iniciales, no se perdió ni se duplicó ninguna partícula.

**El segundo momento es necesario**: `Σ id` sola no distingue el caso en que una partícula fuera reemplazada por otra de identidad complementaria (perder la 3 y duplicar la 5 mantiene la suma si simultáneamente se pierde la 7 y se duplica la 5). Con `Σ id²` esa compensación deja de ser posible en la práctica.

## API

| Función | Rol |
|---|---|
| `migrate_particles` | Reubica las partículas fuera de rango; retorna el nuevo `n_local` |
| `replicate_particles` | `Allgatherv` global; devuelve el `offset` del tramo local |
| `particle_checksum` | `Σ id` y `Σ id²` globales |

## Contrato de `migrate_particles`

- `*p` **puede ser reasignado**: la función libera el buffer viejo y devuelve uno nuevo. `*capacity` se actualiza. No conservar punteros al arreglo a través de la llamada.
- Al retornar, el arreglo local queda **ordenado por clave creciente**.
- El valor de retorno puede ser `0`.
- Colectiva.

## Casos borde

| Caso | Comportamiento |
|---|---|
| `n_local == 0` al entrar | Legal; envía 0 a todos |
| `n_recv == 0` al salir | Legal; se reserva un buffer de 1 elemento para no devolver `NULL` |
| Ninguna partícula migra | `Alltoallv` con todo en la diagonal; correcto y barato |
| Partícula fuera del box congelado | `morton_encode` la acota al borde; migra al proceso del borde |
