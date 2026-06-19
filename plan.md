# Plan de Implementación

## Simulación Paralela del Problema de N-Cuerpos en 3D mediante MPI + OpenMP con Barnes–Hut y Balance Dinámico

---

# 1. Análisis de la arquitectura propuesta

## 1.1 Arquitectura general

La aplicación se implementará como una simulación híbrida MPI/OpenMP:

- **MPI:** distribución de partículas entre procesos.
- **OpenMP:** paralelización intra-proceso.
- **Barnes-Hut:** reducción de complejidad del cálculo de fuerzas.
- **Morton:** particionamiento espacial.
- **LET:** intercambio de información remota.
- **Leapfrog:** integración temporal.

La simulación seguirá el siguiente flujo:

```
Inicialización
    ↓
Generación claves Morton
    ↓
Particionamiento inicial
    ↓
Distribución MPI
    ↓
Bucle temporal:

    Construcción octree local
    ↓
    Construcción LET
    ↓
    Comunicación MPI
    ↓
    Integración LET
    ↓
    Cálculo de fuerzas
    ↓
    Leapfrog
    ↓
    Migración de partículas
    ↓
    Evaluación del balance
    ↓
    Rebalanceo (si corresponde)
```

---

## 1.2 Componentes principales

| Componente           | Función                  |
| -------------------- | ------------------------ |
| Particle System      | Almacenar partículas     |
| Morton Encoder       | Ordenamiento espacial    |
| Domain Decomposition | Particionamiento MPI     |
| Octree Builder       | Construcción del árbol   |
| Barnes-Hut Solver    | Cálculo de fuerzas       |
| LET Manager          | Información remota       |
| MPI Communication    | Intercambio de datos     |
| Migration System     | Movimiento de partículas |
| Load Balancer        | Reparto dinámico         |
| Leapfrog Integrator  | Evolución temporal       |
| Metrics Module       | Medición                 |
| Validation Suite     | Pruebas                  |

---

## 1.3 Dependencias

```
Partículas
    ↓
Morton
    ↓
Particiones
    ↓
Octree
    ↓
LET
    ↓
Barnes-Hut
    ↓
Leapfrog
    ↓
Migración
    ↓
Balance
```

---

## 1.4 Partes técnicamente complejas

### Alta complejidad

- Generación de LET.
- Comunicación distribuida.
- Rebalanceo dinámico.
- Integración árbol local + LET.

### Complejidad media

- Octree.
- Barnes-Hut.
- Migración.

### Baja complejidad

- Leapfrog.
- Morton.
- Instrumentación.

---

# 2. Descomposición en módulos

---

## Módulo: Particle

### Objetivo

Representar las partículas.

### Entrada

Archivo inicial.

### Salida

Vector de partículas.

### Complejidad

Baja.

---

## Módulo: Leapfrog

Responsable de:

- Drift.
- Kick.

Complejidad: Baja.

---

## Módulo: Morton

Generar claves espaciales.

Entrada:

- Posición.

Salida:

- uint64 Morton.

Complejidad: Baja.

---

## Módulo: Domain Partition

Responsable:

- Dividir rango Morton.

Salida:

- Límites por proceso.

Complejidad: Media.

---

## Módulo: Octree

Responsabilidades:

- Insertar partículas.
- Subdividir nodos.
- Calcular centros de masa.

Complejidad: Alta.

---

## Módulo: Barnes-Hut

Responsabilidades:

- Recorrer árbol.
- Aplicar criterio θ.

Complejidad: Alta.

---

## Módulo: LET

Responsabilidades:

- Identificar nodos necesarios.
- Empaquetar.
- Enviar.

Complejidad: Muy alta.

---

## Módulo: MPI

Responsabilidades:

- Allgather.
- Sendrecv.
- Alltoallv.

Complejidad: Media.

---

## Módulo: Migración

Responsabilidades:

- Detectar partículas fuera del dominio.
- Enviar.

Complejidad: Media.

---

## Módulo: Balance

Responsabilidades:

- Medir carga.
- Recalcular límites.

Complejidad: Alta.

---

## Módulo: OpenMP

Responsabilidades:

- Paralelizar fuerzas.

Complejidad: Baja.

---

## Módulo: Métricas

Responsabilidades:

- Timers.
- Estadísticas.

Complejidad: Baja.

---

# 3. Orden recomendado de implementación

---

## Etapa 1

### Objetivo

Versión secuencial.

Implementar:

- Partículas.
- Leapfrog.
- Fuerza O(N²).

Validación:

- Órbitas simples.

---

## Etapa 2

Agregar Barnes-Hut local.

Validar:

- Error respecto O(N²).

---

## Etapa 3

Agregar Morton.

Validar:

- Localidad espacial.

---

## Etapa 4

Agregar octree.

Validar:

- Masa total.
- Centros de masa.

---

## Etapa 5

OpenMP.

Validar:

- Mismo resultado.
- Speedup.

---

## Etapa 6

MPI estático.

Validar:

- Conservación de partículas.

---

## Etapa 7

Migración.

Validar:

- Ninguna partícula perdida.

---

## Etapa 8

LET.

Validar:

- Comparación contra árbol global.

---

## Etapa 9

Rebalanceo.

Validar:

- Reducción del desbalance.

---

# 4. Estructuras de datos

---

## Particle

```cpp
struct Particle {
    double x,y,z;
    double vx,vy,vz;
    double ax,ay,az;
    double mass;
    uint64_t morton;
};
```

---

## OctreeNode

```cpp
struct Node {
    BoundingBox box;

    double mass;

    double cx,cy,cz;

    bool leaf;

    int children[8];

    vector<int> particles;
};
```

---

## LETNode

```cpp
struct LETNode {
    double mass;

    double cx,cy,cz;

    BoundingBox box;

    bool leaf;
};
```

---

## MPIParticle

```cpp
struct MPIParticle {
    double pos[3];
    double vel[3];
    double mass;
};
```

---

## Partition

```cpp
struct Partition {
    uint64_t min_key;
    uint64_t max_key;
};
```

---

## Metrics

```cpp
struct Metrics {
    double compute_time;
    double mpi_time;
    double rebalance_time;
    int particles;
};
```

---

# 5. Diseño del algoritmo distribuido

---

# Construcción del árbol

Cada proceso posee:

- Partículas locales.
- Límites Morton.

Construye únicamente su octree.

---

# Generación del LET

Para cada proceso remoto:

1. Recorrer árbol local.
2. Evaluar criterio:

```
s / d < θ
```

Si se cumple:

- enviar nodo.

Si no:

- descender.

---

# Datos enviados

Cada nodo:

- masa.
- centro de masa.
- tamaño.
- bounding box.

No se envían partículas.

---

# Comunicación MPI

## Recomendadas

- MPI_Allgather.
- MPI_Alltoallv.
- MPI_Isend.
- MPI_Irecv.

---

## Fases

### LET

Cada iteración.

### Migración

Cada iteración.

### Rebalanceo

Cada K iteraciones.

---

# Cálculo de fuerzas

Árbol utilizado:

```
Árbol local
+
LET recibido
```

Para cada partícula:

```
recorrer(arbol_local)
recorrer(LET)
```

Aplicando:

```
s/d < θ
```

---

# Migración

Después del Leapfrog:

```cpp
if particle.key outside local range:
    migrate
```

---

# Rebalanceo

Cada:

```text
20-50 iteraciones
```

Calcular:

```
imbalance =
max(load)/avg(load)
```

Si:

```
imbalance > 1.2
```

Repartir nuevamente.

---

# 6. Estrategia OpenMP

---

## Paralelizar

```cpp
#pragma omp parallel for schedule(dynamic)
for(i=0;i<N;i++)
```

---

## Shared

- Árbol.
- LET.

---

## Private

- Fuerza local.
- Pilas de recorrido.

---

## Riesgos

### Acumulación compartida

Evitar:

```cpp
force +=
```

Cada hilo trabaja sobre una partícula.

---

## Justificación dynamic

Distintas partículas recorren distinta cantidad de nodos.

El costo no es uniforme.

---

# 7. Estrategia de validación

---

## Morton

Entrada:

- Puntos conocidos.

Esperado:

- Orden correcto.

---

## Octree

Verificar:

- Masa total.
- Centro de masa.

---

## Barnes-Hut

Comparar:

- Error relativo.

---

## LET

Comparar:

- Fuerzas obtenidas.
- Árbol global.

---

## MPI

Verificar:

- Suma de partículas.

---

## Migración

Verificar:

- Ninguna partícula perdida.

---

## Balance

Verificar:

- Desbalance disminuye.

---

# 8. Evaluación experimental

---

## Tamaños

| N         |
| --------- |
| 10.000    |
| 50.000    |
| 100.000   |
| 500.000   |
| 1.000.000 |

---

## MPI

- 1
- 2
- 4
- 8
- 16

---

## OpenMP

- 1
- 2
- 4
- 8
- 16

---

## Métricas

- Tiempo total.
- Speedup.
- Eficiencia.
- Tiempo MPI.
- Tiempo LET.
- Tiempo rebalanceo.

---

## Escalabilidad fuerte

N fijo.

---

## Escalabilidad débil

N/P constante.

---

# 9. Cronograma (1 mes, 2 estudiantes)

---

## Semana 1

### Estudiante A

- Partículas.
- Leapfrog.
- Morton.

### Estudiante B

- Fuerza O(N²).
- Validación.

### Entregable

Versión secuencial.

---

## Semana 2

A:

- Octree.

B:

- Barnes-Hut.

Entregable:

Barnes-Hut secuencial.

---

## Semana 3

A:

- MPI.
- Particiones.

B:

- OpenMP.
- Migración.

Entregable:

Versión híbrida básica.

---

## Semana 4

A:

- LET.

B:

- Balance.

Ambos:

- Experimentos.
- Informe.

Entregable:

Sistema final.

---

# 10. Riesgos técnicos

| Riesgo              | Impacto | Probabilidad | Mitigación                 |
| ------------------- | ------- | ------------ | -------------------------- |
| LET incorrecto      | Alto    | Alta         | Comparar con árbol global  |
| Partículas perdidas | Alto    | Media        | Verificaciones globales    |
| Desbalance excesivo | Medio   | Alta         | Rebalanceo                 |
| Deadlocks MPI       | Alto    | Media        | Comunicación no bloqueante |
| Error Barnes-Hut    | Medio   | Media        | Comparación O(N²)          |
| Sobrecarga LET      | Medio   | Alta         | Instrumentación            |

---

# 11. Producto mínimo viable (MVP)

El MVP recomendado es:

- MPI estático.
- Morton.
- Octree local.
- Barnes-Hut.
- OpenMP.
- Migración.
- Sin rebalanceo.
- LET simplificado.

Configuración:

- 2–4 procesos MPI.
- 4–8 hilos OpenMP.
- 100.000 partículas.

Este sistema ya permite medir:

- Speedup.
- Escalabilidad.
- Comunicación.
- Beneficio de Barnes-Hut.

---

## Funcionalidades postergables

Si el tiempo resulta insuficiente:

1. Rebalanceo dinámico.
2. LET optimizado.
3. Comunicación asíncrona.
4. Escalabilidad fuerte de gran tamaño.
5. Métricas avanzadas.

El orden recomendado de prioridades es:

1. Barnes-Hut local.
2. MPI estático.
3. OpenMP.
4. Migración.
5. LET.
6. Balance dinámico.

De esta forma se garantiza obtener resultados experimentales válidos incluso si las funcionalidades más complejas no llegan a completarse.
