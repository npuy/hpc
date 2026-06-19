# Plan de Implementacion

## Simulacion Paralela del Problema de N-Cuerpos en 3D mediante MPI + OpenMP con Barnes-Hut y Balance Dinamico

## Arquitectura general

Simulacion hibrida MPI/OpenMP:

- **MPI:** distribucion de particulas entre procesos.
- **OpenMP:** paralelizacion intra-proceso.
- **Barnes-Hut:** reduccion de complejidad del calculo de fuerzas.
- **Morton:** particionamiento espacial.
- **LET:** intercambio de informacion remota.
- **Leapfrog:** integracion temporal.

### Flujo de la simulacion

```
Inicializacion
    |
Generacion claves Morton
    |
Particionamiento inicial
    |
Distribucion MPI
    |
Bucle temporal:
    Construccion octree local
    -> Construccion LET
    -> Comunicacion MPI
    -> Integracion LET
    -> Calculo de fuerzas
    -> Leapfrog
    -> Migracion de particulas
    -> Evaluacion del balance
    -> Rebalanceo (si corresponde)
```

## Componentes principales

| Componente           | Funcion                  |
| -------------------- | ------------------------ |
| Particle System      | Almacenar particulas     |
| Morton Encoder       | Ordenamiento espacial    |
| Domain Decomposition | Particionamiento MPI     |
| Octree Builder       | Construccion del arbol   |
| Barnes-Hut Solver    | Calculo de fuerzas       |
| LET Manager          | Informacion remota       |
| MPI Communication    | Intercambio de datos     |
| Migration System     | Movimiento de particulas |
| Load Balancer        | Reparto dinamico         |
| Leapfrog Integrator  | Evolucion temporal       |
| Metrics Module       | Medicion                 |
| Validation Suite     | Pruebas                  |

## Cadena de dependencias

```
Particulas -> Morton -> Particiones -> Octree -> LET -> Barnes-Hut -> Leapfrog -> Migracion -> Balance
```

## Complejidad por componente

- **Alta:** LET, comunicacion distribuida, rebalanceo dinamico, integracion arbol local + LET
- **Media:** Octree, Barnes-Hut, migracion
- **Baja:** Leapfrog, Morton, instrumentacion

## Estructuras de datos

### Particle

```c
typedef struct {
    double pos[3];      // x, y, z
    double vel[3];      // vx, vy, vz
    double acc[3];      // ax, ay, az
    double mass;
    uint64_t morton;
} Particle;
```

### OctreeNode

```c
typedef struct {
    double min[3], max[3]; // bounding box
    double mass;
    double cm[3];          // centro de masa
    int is_leaf;
    int children[8];
    int *particles;
    int n_particles;
} OctreeNode;
```

### LETNode

```c
typedef struct {
    double mass;
    double cm[3];
    double min[3], max[3];
    int is_leaf;
} LETNode;
```

### Partition

```c
typedef struct {
    uint64_t min_key;
    uint64_t max_key;
} Partition;
```

### Metrics

```c
typedef struct {
    double compute_time;
    double mpi_time;
    double rebalance_time;
    int particles;
} Metrics;
```

## Etapas de implementacion

| Etapa | Objetivo                | Validacion                        |
| ----- | ----------------------- | --------------------------------- |
| 1     | Version secuencial      | Orbitas simples                   |
| 2     | Barnes-Hut local        | Error respecto O(N^2)             |
| 3     | Morton                  | Localidad espacial                |
| 4     | Octree                  | Masa total, centros de masa       |
| 5     | OpenMP                  | Mismo resultado, speedup          |
| 6     | MPI estatico            | Conservacion de particulas        |
| 7     | Migracion               | Ninguna particula perdida         |
| 8     | LET                     | Comparacion contra arbol global   |
| 9     | Rebalanceo              | Reduccion del desbalance          |

## Cronograma (1 mes, 2 estudiantes)

### Semana 1 — Version secuencial

- **A:** Particulas, Leapfrog, Morton
- **B:** Fuerza O(N^2), Validacion

### Semana 2 — Barnes-Hut secuencial

- **A:** Octree
- **B:** Barnes-Hut

### Semana 3 — Version hibrida basica

- **A:** MPI, Particiones
- **B:** OpenMP, Migracion

### Semana 4 — Sistema final

- **A:** LET
- **B:** Balance
- **Ambos:** Experimentos, Informe

## Evaluacion experimental

- **Tamanos:** 10K, 50K, 100K, 500K, 1M particulas
- **MPI:** 1, 2, 4, 8, 16 procesos
- **OpenMP:** 1, 2, 4, 8, 16 hilos
- **Metricas:** Tiempo total, speedup, eficiencia, tiempo MPI/LET/rebalanceo
- **Escalabilidad fuerte:** N fijo, variar procesos
- **Escalabilidad debil:** N/P constante

## MVP (Producto Minimo Viable)

- MPI estatico + Morton + Octree local + Barnes-Hut + OpenMP + Migracion
- Sin rebalanceo, LET simplificado
- Config: 2-4 MPI, 4-8 OpenMP, 100K particulas

## Prioridades si falta tiempo

1. Barnes-Hut local
2. MPI estatico
3. OpenMP
4. Migracion
5. LET
6. Balance dinamico

## Riesgos tecnicos

| Riesgo              | Impacto | Probabilidad | Mitigacion                 |
| ------------------- | ------- | ------------ | -------------------------- |
| LET incorrecto      | Alto    | Alta         | Comparar con arbol global  |
| Particulas perdidas | Alto    | Media        | Verificaciones globales    |
| Desbalance excesivo | Medio   | Alta         | Rebalanceo                 |
| Deadlocks MPI       | Alto    | Media        | Comunicacion no bloqueante |
| Error Barnes-Hut    | Medio   | Media        | Comparacion O(N^2)         |
| Sobrecarga LET      | Medio   | Alta         | Instrumentacion            |
