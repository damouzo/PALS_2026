# Comparativa de Workflow Managers en Bioinformática

> *Navegar el ecosistema: ¿qué herramienta elegir y por qué?*

---

## ¿Por qué existen tantos workflow managers?

Los workflow managers en bioinformática no surgieron de un diseño centralizado, sino de necesidades independientes en diferentes comunidades. Todos resuelven el mismo problema fundamental — **gestionar pipelines de análisis complejos, reproducibles y escalables** — pero con filosofías y syntaxis muy distintas.

Entender el ecosistema es importante porque:
1. Colaborarás con personas que usan otras herramientas.
2. Encontrarás pipelines publicados en distintos lenguajes.
3. Podrás tomar decisiones informadas sobre qué adoptar en tu grupo.

---

## El ecosistema global de workflow managers

### Principales herramientas en bioinformática

| Herramienta | Año | Lenguaje base | Paradigma | Popularidad |
|-------------|-----|---------------|-----------|-------------|
| **Snakemake** | 2012 | Python | Basado en reglas (Make) | ⭐⭐⭐⭐⭐ (academia) |
| **Nextflow** | 2013 | Groovy/DSL propio | Flujo de datos (dataflow) | ⭐⭐⭐⭐⭐ (industria + academia) |
| **WDL** | 2014 | Propio (Broad) | Basado en tareas | ⭐⭐⭐ (Broad/Terra) |
| **CWL** | 2014 | YAML/JSON | Estándar abierto | ⭐⭐ (estándar) |
| **Luigi** | 2012 | Python | Basado en tareas | ⭐⭐ (legado) |
| **Galaxy** | 2005 | Python | GUI + scripting | ⭐⭐⭐ (no-code) |

### Herramientas de propósito general (a veces usadas en bio)

- **Apache Airflow** — orquestación de datos empresariales, también usado en bioinformática
- **Prefect / Dagster** — modernos, orientados a datos, menos adopción en bio
- **Cromwell** — motor de ejecución para WDL, desarrollado por el Broad Institute

---

## Snakemake en profundidad

### Origen y filosofía

Snakemake fue creado por **Johannes Köster** en 2012 como parte de su doctorado. El nombre ya dice todo: es un homenaje a **GNU Make**, el sistema de construcción de Unix, pero con la potencia de Python.

La idea central es **pensar en los outputs en lugar de los inputs**: "¿qué fichero quiero generar? ¿qué necesito para generarlo?". El motor de Snakemake luego calcula automáticamente el DAG hacia atrás.

### Sintaxis y conceptos clave

```python
# Snakefile — Pipeline RNA-seq simple en Snakemake

# Configuración
configfile: "config/config.yaml"

SAMPLES = config["samples"]

# Regla final (objetivo): lo que queremos generar
rule all:
    input:
        expand("results/counts/{sample}_counts.txt", sample=SAMPLES),
        "results/multiqc/multiqc_report.html"

# Regla de QC
rule fastqc:
    input:
        r1 = "data/fastq/{sample}_R1.fastq.gz",
        r2 = "data/fastq/{sample}_R2.fastq.gz"
    output:
        html_r1 = "results/fastqc/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc/{sample}_R2_fastqc.html"
    threads: 4
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -t {threads} {input.r1} {input.r2} -o results/fastqc/"

# Regla de alineamiento
rule star_align:
    input:
        r1 = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        genome = config["genome_index"]
    output:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.genome} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix results/star/{wildcards.sample}/ \
             --outSAMtype BAM SortedByCoordinate
        """
```

### Características destacadas de Snakemake

- **Wildcards**: el equivalente a las variables de samples. `{sample}` en los paths se infiere automáticamente.
- **Integración con conda**: `conda:` en cada regla especifica el entorno. Automáticamente creado con `--use-conda`.
- **Re-ejecución inteligente**: solo re-ejecuta reglas cuyo output no existe o es más antiguo que el input.
- **Basado en Python**: puedes usar código Python arbitrario en el Snakefile.
- **`snakemake-wrappers`**: repositorio comunitario de reglas pre-escritas, similar a los módulos de nf-core.

### Comunidad y ecosistema

- **Snakemake Workflow Catalog**: repositorio de workflows curados ([snakemake.github.io/snakemake-workflow-catalog](https://snakemake.github.io/snakemake-workflow-catalog/))
- Muy popular en grupos de bioinformática europeos
- Excelente integración con entornos académicos donde Python es el lenguaje dominante

---

## Nextflow en profundidad

### Origen y filosofía

Nextflow fue creado por **Paolo Di Tommaso** en 2013 en el CRG (Centre for Genomic Regulation) de Barcelona. A diferencia de Snakemake, Nextflow adopta el modelo de **programación de flujo de datos (dataflow)**: en lugar de pensar en ficheros de entrada/salida, piensas en **canales** por los que fluyen los datos entre procesos.

El lenguaje se basa en **Groovy** (JVM), aunque en la práctica se usa el DSL propio de Nextflow que abstrae la mayor parte de Groovy. Desde **DSL2** (2020), la sintaxis es modular y componible.

### Sintaxis y conceptos clave

```nextflow
// main.nf — Pipeline RNA-seq en Nextflow (DSL2)

nextflow.enable.dsl=2

params.samples_tsv  = "config/samples.tsv"
params.genome_index = "data/genome/star_index"
params.outdir       = "results"

// Definición de proceso
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    conda "bioconda::fastqc=0.12.1"
    // o container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}")

    script:
    """
    fastqc -t 4 ${r1} ${r2}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}*.bam")

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${genome_index} \
         --readFilesIn ${r1} ${r2} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_ \
         --outSAMtype BAM SortedByCoordinate
    """
}

// Workflow
workflow {
    // Crear canal desde TSV
    ch_samples = Channel
        .fromPath(params.samples_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.fastq_R1), file(row.fastq_R2)) }

    // Ejecutar procesos
    FASTQC(ch_samples)
    STAR_ALIGN(ch_samples, params.genome_index)
}
```

### Características destacadas de Nextflow

- **Canales (Channels)**: el eje central del paradigma. Los datos "fluyen" por los procesos.
- **`-resume`**: re-ejecuta solo los procesos que han cambiado, usando un cache de tareas.
- **Portabilidad sin cambios de código**: el mismo `main.nf` corre en local, SLURM, AWS Batch, Google Cloud, Azure.
- **`nextflow.config`**: configuración separada del código para adaptar a distintos entornos.
- **Contenedores nativos**: Docker, Singularity, Podman integrados directamente.
- **nf-core**: ecosistema de más de 100 pipelines de producción, bien documentados y mantenidos.

### nf-core: el ecosistema de pipelines curados

nf-core es un proyecto comunitario que:
- Mantiene pipelines de "producción" para los análisis más comunes (RNA-seq, scRNA-seq, ATAC-seq, variant calling, etc.)
- Establece un estándar de buenas prácticas para desarrollar en Nextflow
- Ofrece módulos reutilizables para las herramientas bioinformáticas más comunes
- Con `nf-core/rnaseq` puedes correr un análisis completo de RNA-seq con una sola línea

```bash
# Ejecutar el pipeline nf-core/rnaseq
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile singularity
```

---

## Snakemake vs. Nextflow: la comparativa central

### Filosofía

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **Paradigma** | Basado en reglas (outputs → inputs) | Basado en flujo de datos (canales) |
| **Inspiración** | GNU Make + Python | Erlang + Groovy |
| **Forma de pensar** | "¿Qué ficheros quiero crear?" | "¿Cómo fluyen los datos entre procesos?" |
| **Curva de aprendizaje** | Suave (si sabes Python) | Media (nuevo DSL, concepto de canales) |

### Lenguaje y código

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **Lenguaje base** | Python | Groovy (JVM) |
| **DSL propio** | Mínimo (extiende Python) | Sí (DSL2 propio) |
| **Integración con Python** | Total (código Python en Snakefile) | Limitada (Nextflow no es Python) |
| **Integración con R** | Vía shell | Vía shell |
| **Debugging** | Python-native (excepciones claras) | Más complejo (logs de tarea) |

### Ejecución y portabilidad

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **HPC (SLURM, SGE)** | Sí (profiles) | Sí (executors) |
| **Cloud** | Sí (AWS, GCP, Azure) | Sí (AWS Batch, GCP, Azure) |
| **Kubernetes** | Sí | Sí |
| **Local** | Sí | Sí |
| **Portabilidad entre entornos** | Buena (profiles) | Excelente (nextflow.config) |

### Contenedores y entornos

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **Conda** | Integrado (por regla) | Integrado (por proceso) |
| **Docker** | Sí | Sí |
| **Singularity** | Sí | Sí (más maduro) |
| **Entornos por herramienta** | Sí (cada regla tiene su env) | Sí (cada proceso tiene su container) |

### Re-ejecución y caché

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **Mecanismo** | Timestamps de ficheros | Hash de inputs + código |
| **Cache explícito** | No | Sí (`.nextflow/cache`) |
| **Re-run selectivo** | Basado en ficheros modificados | `-resume` (hash-based) |
| **Robustez** | Buena | Excelente |

### Comunidad y ecosistema

| Aspecto | Snakemake | Nextflow |
|---------|-----------|----------|
| **Comunidad principal** | Academia europea | Academia + industria global |
| **Pipelines curados** | Snakemake Workflow Catalog | nf-core (más de 100 pipelines) |
| **Publicaciones** | ~7,000 citaciones | ~7,000+ citaciones |
| **Empresa detrás** | Comunidad open source | Seqera Labs (comercial) |
| **Soporte comercial** | No | Sí (Seqera Platform) |

### ¿Cuándo elegir uno u otro?

**Elige Snakemake si:**
- Tu grupo está orientado a Python y quieres la menor fricción.
- Tus análisis son explorativos y los outputs son ficheros bien definidos.
- Prefieres un paradigma más cercano a Make (declarativo, orientado a ficheros).
- Trabajas principalmente en entornos académicos europeos donde Snakemake es la norma.

**Elige Nextflow si:**
- Quieres acceder al ecosistema nf-core y aprovechar pipelines de producción ya hechos.
- Tus pipelines procesan flujos de datos complejos (muchos samples, datos heterogéneos).
- Necesitas portabilidad real entre HPC, cloud y local sin tocar el código.
- Tu grupo o institución tiene un perfil más "industria" o trabajas con equipos grandes.
- Quieres re-ejecuciones con `-resume` verdaderamente robustas (hash-based).

**La realidad en 2024-2026:**
En la práctica, ambas herramientas son excelentes y la elección a menudo depende más de quién está en tu grupo, qué usa tu institución y qué pipelines necesitas, que de diferencias técnicas fundamentales.

---

## WDL (Workflow Description Language)

Desarrollado por el **Broad Institute** y usado principalmente en **Terra** (plataforma de análisis del Broad/Google).

```wdl
# Ejemplo WDL
version 1.0

task FastQC {
    input {
        File fastq_r1
        File fastq_r2
    }
    command <<<
        fastqc ~{fastq_r1} ~{fastq_r2}
    >>>
    output {
        Array[File] html_reports = glob("*_fastqc.html")
    }
    runtime {
        docker: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
        memory: "4 GB"
        cpu: 4
    }
}
```

**Contexto**: WDL es la elección natural si trabajas con datos de TCGA, GNOMAD o con colaboradores del Broad/NIH. Fuera de ese ecosistema, tiene menos adopción.

---

## CWL (Common Workflow Language)

Un estándar open source diseñado para la interoperabilidad entre plataformas.

**Ventaja**: Los pipelines CWL pueden ejecutarse en múltiples motores (cwltool, Arvados, Galaxy, etc.).  
**Desventaja**: Extremadamente verboso y difícil de escribir manualmente.

**Contexto**: Usado principalmente en consorcios y proyectos que requieren interoperabilidad estricta. Poco práctico para el día a día.

---

## Galaxy

Plataforma web con interfaz gráfica orientada a usuarios sin programación.

**Ventaja**: Accesible para clínicos, biólogos sin programación, educación.  
**Desventaja**: Difícil de versionar, limita la personalización, no es adecuado para pipelines de producción a gran escala.

**Contexto**: Excelente para demostraciones y para colaboradores no bioinformáticos. Tiene una comunidad activa y muchas herramientas instaladas.

---

## Resumen ejecutivo

```
Accesibilidad (no-código)
         ↑
    Galaxy
         |
    Snakemake ← (buen punto de entrada si ya sabes Python)
         |
    Nextflow  ← (mayor potencial de escala, ecosistema nf-core)
         |
    WDL/CWL
         ↓
Flexibilidad y escala
```

> **Para el contexto de un summer school de bioinformática:** Nextflow y Snakemake son las opciones más relevantes. Ambas son herramientas del estado del arte. La diferencia práctica más importante para un predoc es el ecosistema: **nf-core** proporciona acceso inmediato a pipelines de producción que pueden acelerar enormemente el análisis exploratorio antes de construir pipelines propios.

---

## Referencias

- Köster, J. & Rahmann, S. *Snakemake — a scalable bioinformatics workflow engine.* Bioinformatics 28, 2520–2522 (2012)
- Di Tommaso, P. et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology 35, 316–319 (2017)
- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core](https://nf-co.re/)
- [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
