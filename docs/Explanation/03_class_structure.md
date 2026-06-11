# Estructura del Seminario: Workflow Management en Bioinformática

> *Plan docente para el Summer School de Bioinformática*  
> **Audiencia**: Predocs y postdocs tempranos  
> **Duración estimada**: 3–4 horas (sesión completa + hands-on)  
> **Herramienta principal**: Nextflow

---

## Objetivo general del seminario

Al finalizar esta sesión, los asistentes entenderán por qué un workflow manager es una herramienta esencial para el bioinformático moderno, serán capaces de identificar los conceptos clave de Nextflow, y habrán dado sus primeros pasos prácticos desarrollando un pipeline básico.

---

## Perfil de la audiencia

**Quiénes son:**
- Estudiantes de doctorado en etapas tempranas-medias y postdocs iniciando
- Biologías moleculares, bioinformáticos, clínicos computacionales
- Rango amplio de experiencia en programación: desde "sé lo básico de R" hasta "hago Python a diario"
- Probable experiencia con análisis de RNA-seq o scRNA-seq, al menos como usuarios de pipelines o scripts
- Frustrados con análisis no reproducibles, scripts caóticos o dependencias de un único cluster

**Qué necesitan:**
- Argumentos para convencer a su supervisor de adoptar nuevas herramientas
- Herramientas concretas y no solo teoría
- Confianza para empezar a usar estas tecnologías por su cuenta
- Recursos a los que volver después del curso

---

## Estructura del seminario

### Bloque 1 — El problema: ¿por qué necesitamos workflow management?
**Duración estimada:** 30–40 minutos  
**Formato:** Presentación + discusión participativa

#### Capítulo 1.1 — La historia que todos conocemos
- Empezar con una pregunta directa al aula: *"¿Quién alguna vez no ha podido reproducir un análisis propio de hace 6 meses?"*
- Recorrido narrativo por la evolución de los workflows (desde la terminal hasta el workflow manager)
- Mostrar la progresión de complejidad vs. organización de forma visual
- **Mensaje clave:** No hay un salto abrupto: es una evolución natural y todos estáis ya en algún punto de ese camino

#### Capítulo 1.2 — ¿Qué es reproducibilidad y por qué importa?
- Definición práctica: otro investigador puede obtener los mismos resultados usando los mismos datos
- El problema de la crisis de reproducibilidad en ciencias de la vida
- La diferencia entre reproducibilidad *técnica* (mismo pipeline, mismos datos → mismos resultados) y *biológica*
- Ejemplos reales de consecuencias: retractaciones, revisiones fallidas, tiempo perdido
- **Mensaje clave:** La reproducibilidad no es solo una virtud académica: es una herramienta de productividad personal

#### Capítulo 1.3 — Los problemas concretos del bioinformático en 2024
- "Funciona en mi máquina" — el problema del entorno
- Versiones de software: DESeq2 2.1 vs DESeq2 1.4 pueden dar resultados distintos
- El problema de escalar: de 5 muestras a 500
- El problema de los colaboradores: enviar el pipeline a otro grupo
- **Mensaje clave:** Estos problemas tienen solución estándar, no tienes que reinventar la rueda

---

### Bloque 2 — El ecosistema: panorama de workflow managers
**Duración estimada:** 20–25 minutos  
**Formato:** Presentación comparativa

#### Capítulo 2.1 — ¿Qué es un workflow manager?
- Definición: software que orquesta, paraleliza, y gestiona la ejecución de pipelines multi-step
- Componentes universales: procesos/reglas, inputs/outputs, dependencias, ejecución distribuida
- El concepto de DAG (Directed Acyclic Graph): visualizar el análisis como un grafo
- **Mensaje clave:** Un workflow manager no escribe el análisis por ti, pero lo ejecuta de forma inteligente

#### Capítulo 2.2 — El mapa del ecosistema
- Rápida panorámica de las opciones: Snakemake, Nextflow, WDL, CWL, Galaxy
- Criterios para elegir: comunidad, ecosistema de pipelines, curva de aprendizaje, entorno de trabajo
- Posicionar Snakemake vs. Nextflow como las dos opciones dominantes en academia 2024-2026
- Honestidad sobre las limitaciones: ninguna herramienta es perfecta para todos los casos
- **Mensaje clave:** Nextflow tiene el ecosistema más rico para bioinformática en este momento gracias a nf-core

#### Capítulo 2.3 — nf-core: por qué es un game-changer
- ¿Qué es nf-core? Community-driven pipelines de producción para bioinformática
- Pipelines más usados: `rnaseq`, `sarek` (variant calling), `scrnaseq`, `atacseq`, `chipseq`
- La promesa: `nextflow run nf-core/rnaseq` —uno de los análisis más reproducibles disponibles
- El ecosistema de módulos: reutilizar en lugar de re-escribir
- **Mensaje clave:** nf-core es el punto de acceso más directo al valor de Nextflow para un bioinformático

---

### Bloque 3 — Nextflow: conceptos fundamentales
**Duración estimada:** 40–50 minutos  
**Formato:** Presentación + live coding

#### Capítulo 3.1 — Arquitectura y filosofía de Nextflow
- El modelo de dataflow: los datos "fluyen" entre procesos a través de canales
- Comparación visual: el pipeline de bash secuencial vs. el DAG de Nextflow
- Los tres grandes elementos: **Processes**, **Channels**, **Workflows**
- La separación entre lógica de análisis (`main.nf`) y configuración de ejecución (`nextflow.config`)
- **Mensaje clave:** Un proceso de Nextflow describe *qué hace* una herramienta; el canal describe *con qué datos* lo hace

#### Capítulo 3.2 — Anatomía de un proceso
- `input:`, `output:`, `script:`/`exec:` — los bloques esenciales
- `tag`, `publishDir`, `label` — las directivas más importantes
- El concepto de task: cada proceso × cada sample = una task independiente
- El directorio `work/`: cómo Nextflow aísla cada ejecución
- El hash hexadecimal: identificador único de cada task
- **Mensaje clave:** Cada task es un proceso aislado con su propio entorno y archivos

#### Capítulo 3.3 — Channels: el corazón del dataflow
- `Channel.of()`, `Channel.fromPath()`, `Channel.fromFilePairs()`
- Tuples para combinar datos (sample ID + archivos)
- Operators básicos: `.map()`, `.filter()`, `.collect()`, `.flatten()`
- El concepto de metadata map: la forma idiomática de manejar sample info en nf-core
- **Mensaje clave:** Los canales son la clave para la paralelización automática

#### Capítulo 3.4 — Resume, cache y reproducibilidad
- La magia de `-resume`: por qué es diferente a otros sistemas
- El hash de tasks: cómo Nextflow decide qué re-ejecutar
- La importancia del directorio `work/`
- Cuándo *no* usar `-resume` (cambios en el código de un proceso)
- **Mensaje clave:** `-resume` es una de las funcionalidades más valiosas del workflow manager

#### Capítulo 3.5 — Configuración y portabilidad
- `nextflow.config`: el archivo que hace el pipeline portátil
- Profiles: `local`, `slurm`, `hpc_generic`, `singularity`, `docker`
- Ejecutores: local, SLURM, SGE, AWS Batch, Google Cloud
- Contenedores: por qué un proceso con imagen Docker es la máxima reproducibilidad
- **Mensaje clave:** El mismo `main.nf` sin tocar puede correr en tu portátil, el HPC de tu institución y AWS

---

### Bloque 4 — Hands-on: tus primeros pasos con Nextflow
**Duración estimada:** 60–90 minutos  
**Formato:** Workshop interactivo

#### Capítulo 4.1 — Configuración del entorno (15 min)
- Instalación de Nextflow (1 línea, solo requiere Java)
- Opciones de entorno: local, GitHub Codespaces, cluster HPC
- Verificar instalación: `nextflow -version`
- Directorio de trabajo y estructura del proyecto

#### Capítulo 4.2 — Hello World: el primer pipeline (20 min)
- Explicar el pipeline hello world oficial
- Lanzar: `nextflow run main.nf`
- Explorar el directorio `work/`
- Explorar los resultados
- **Ejercicio:** Modificar el parámetro de entrada y ver cómo cambia el comportamiento

#### Capítulo 4.3 — Resume en acción (15 min)
- Modificar el proceso `convert_to_upper` para usar `rev`
- Re-lanzar con `-resume`
- Ver qué tasks se cachean y cuáles se re-ejecutan
- **Ejercicio:** Cambiar el parámetro de entrada y observar el comportamiento del cache

#### Capítulo 4.4 — De Hello World a bioinformática (30 min)
- Adaptar el esquema a un análisis de RNA-seq simplificado
- Channels con sample metadata (tuples: `[sample_id, fastq_r1, fastq_r2]`)
- Proceso FASTQC simple
- `publishDir` para copiar resultados
- **Ejercicio:** Añadir un segundo proceso que depende del primero

#### Capítulo 4.5 — ¿Cómo ejecutar un pipeline nf-core? (15 min)
- `nextflow run nf-core/rnaseq --help`
- El samplesheet CSV de nf-core
- Un run de test: `nextflow run nf-core/rnaseq -profile test,docker`
- Explicar cómo adaptar un run real

---

### Bloque 5 — Más allá del Hello World: scalability y mejores prácticas
**Duración estimada:** 30 minutos  
**Formato:** Presentación + Q&A

#### Capítulo 5.1 — Modules: reutilizando en lugar de re-escribiendo
- El concepto de módulo DSL2: un proceso en su propio archivo
- `include { FASTQC } from './modules/fastqc'`
- Módulos de nf-core: `nf-core modules list`
- La filosofía de "don't repeat yourself" en Nextflow
- **Mensaje clave:** Construir sobre módulos existentes acelera enormemente el desarrollo

#### Capítulo 5.2 — Recursos, labels y configuración HPC
- La directiva `label` para asignar recursos
- `process.cpus`, `process.memory`, `process.time`
- Dynamic resources: reintentar con más memoria si falla
- Configuración de SLURM: el ejecutor más común en academia
- **Mensaje clave:** Nextflow gestiona los recursos automáticamente si le dices cómo

#### Capítulo 5.3 — Control de versiones y buenas prácticas
- Estructura de proyecto recomendada: `main.nf`, `nextflow.config`, `modules/`, `conf/`, `docs/`
- Git para versionar el pipeline (no los datos)
- `nextflow.config` en el repositorio vs. configuraciones locales (ignoradas por git)
- Cómo citar y documentar versiones de software
- **Mensaje clave:** Un pipeline sin Git no es un pipeline de producción

#### Capítulo 5.4 — Recursos para seguir aprendiendo
- [training.nextflow.io](https://training.nextflow.io/) — La mejor fuente oficial
- [nf-core Slack](https://nf-co.re/join/slack) — Comunidad muy activa y acogedora
- [nf-core bytesize talks](https://nf-co.re/events/bytesize) — Charlas cortas sobre temas específicos
- GitHub de pipelines nf-core como referencia de buenas prácticas
- Canal de YouTube de Nextflow

---

## Mensajes clave del seminario

Estos son los 10 mensajes que un asistente debería llevarse al final:

1. **La reproducibilidad no es opcional** — es una herramienta de productividad personal y un estándar de calidad científica.
2. **La evolución hacia los workflow managers es natural** — no es un salto abrupto; ya estás en ese camino.
3. **El DAG es la clave conceptual** — pensar en flujos de datos, no en secuencias de pasos.
4. **Nextflow separa código de configuración** — el mismo pipeline, distintos entornos.
5. **`-resume` cambia todo** — la capacidad de re-ejecutar solo lo que ha cambiado elimina uno de los mayores costos del desarrollo de pipelines.
6. **No construyas desde cero** — nf-core tiene probablemente el análisis que necesitas, o al menos los módulos.
7. **Los contenedores son la máxima reproducibilidad** — si tu proceso tiene una imagen Docker, es reproducible en cualquier máquina con Docker o Singularity.
8. **Git + Nextflow = pipeline de producción** — versionar el código es tan importante como el código mismo.
9. **La curva de aprendizaje vale la inversión** — lo que tardas en aprender Nextflow lo recuperas en el primer análisis a gran escala.
10. **La comunidad nf-core es excepcional** — el Slack de nf-core es uno de los recursos más valiosos en bioinformática computacional en este momento.

---

## Materiales necesarios para el seminario

### Para el instructor
- [ ] Slides de los bloques 1, 2 y 3
- [ ] Entorno con Nextflow instalado para live coding
- [ ] Script hands-on (ver documento `05_nextflow_hands_on.md`)
- [ ] Acceso a internet para `nextflow run nf-core/rnaseq -profile test`

### Para los asistentes
- [ ] Laptop con Java ≥ 11 instalado (o acceso a GitHub Codespaces)
- [ ] Nextflow instalado: `curl -s https://get.nextflow.io | bash`
- [ ] Docker o Singularity (para el ejercicio de nf-core)
- [ ] Cuenta de GitHub (para GitHub Codespaces como alternativa)

---

## Evaluación y feedback

Al final del seminario, recoger:
1. ¿Qué concepto te resultó más difícil de entender?
2. ¿Hay algún aspecto de tu trabajo diario que crees que podrías mejorar con Nextflow?
3. ¿Qué te llevas como primer paso concreto?
4. ¿Qué te gustaría haber visto que no se cubrió?
