# La Evolución de los Workflows en Bioinformática

> *De la terminal al workflow manager: un viaje hacia la reproducibilidad y la escalabilidad*

---

## Introducción

Todo bioinformático ha pasado por las mismas etapas. Empiezas con un par de comandos en la terminal, y antes de que te des cuenta estás gestionando pipelines complejos con decenas de pasos, miles de muestras y análisis que deben ser completamente reproducibles. Este documento traza ese viaje, desde los primeros pasos hasta la adopción de un sistema como Nextflow.

La narrativa sirve un propósito doble: (1) los estudiantes se reconocerán en estas etapas y (2) entenderán *por qué* la complejidad adicional de un workflow manager tiene sentido y es necesaria.

---

## Etapa 0 — "Funciona en mi terminal"

### Descripción

El punto de partida de cualquier bioinformático. Estás siguiendo un tutorial de RNA-seq del paper de STAR, de la documentación de DESeq2 o de un blog post de Bioconductor. Copias los comandos uno a uno en la terminal y ves cómo ocurre la magia.

```bash
# Ejemplo típico: alineamiento con STAR
STAR --runThreadN 8 \
     --genomeDir /data/genome/star_index \
     --readFilesIn sample1_R1.fastq.gz sample1_R2.fastq.gz \
     --outFileNamePrefix results/sample1_ \
     --outSAMtype BAM SortedByCoordinate
```

### Problemas

- **No hay registro** de qué comandos ejecutaste exactamente ni en qué orden.
- Si el servidor se reinicia a mitad del proceso, empiezas de cero.
- Imposible saber exactamente qué versión del software usaste meses después.
- No puedes reproducir el análisis.
- Un error tipográfico puede sobrescribir datos importantes.

### La trampa cognitiva

En esta etapa, la barrera para empezar es mínima. Pero cada hora que ahorras no organizándote se convierte en días perdidos cuando quieres reproducir, compartir o escalar tu análisis.

---

## Etapa 1 — El documento de comandos (el "cuaderno sucio")

### Descripción

Aprendes que deberías guardar los comandos en algún sitio. Creas un archivo de texto, Word, o Notion donde pegas los comandos que vas ejecutando. Sin estructura, sin orden, con partes tachadas y notas al margen.

```
# Notas análisis RNA-seq Proyecto X
# Fecha: algún día de octubre

star comando (el que funcionó al final):
STAR --runThreadN 8 --genomeDir /data/genome/star_index --readFilesIn sample1_R1.fastq.gz...

# OJO: el de arriba no funciona con paired-end, usar este:
STAR --runThreadN 8 --genomeDir /data/genome/star_index --readFilesIn sample1_R1.fastq.gz sample1_R2.fastq.gz...

featureCounts... (buscar parámetros correctos)
```

### Mejoras respecto a la Etapa 0

- Al menos tienes *algo* escrito.
- Puedes buscar comandos que usaste antes.

### Problemas

- El documento es caótico y difícil de seguir.
- Los paths son absolutos y específicos a tu máquina.
- Variables hardcodeadas por todas partes.
- No distingues entre lo que probaste y lo que funcionó.
- Sigue siendo imposible ejecutarlo como pipeline.

---

## Etapa 2 — El primer script real

### Descripción

Alguien (un colega, tu supervisor, un tutorial) te explica que deberías guardar tu código en un **script** con estructura. Creas tu primer `.sh` o `.R` con encabezado, secciones y comentarios.

```bash
#!/bin/bash
# =============================================================
# Script: rnaseq_alignment.sh
# Descripción: Alineamiento de muestras RNA-seq con STAR
# Autor: Tu Nombre
# Fecha: 2024-03-15
# =============================================================

# --- CONFIGURACIÓN ---
THREADS=8
GENOME_DIR="/data/genome/star_index"
INPUT_DIR="/data/fastq"
OUTPUT_DIR="/results/star"

# --- ALINEAMIENTO ---
echo "Iniciando alineamiento de sample1..."
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_DIR \
     --readFilesIn ${INPUT_DIR}/sample1_R1.fastq.gz ${INPUT_DIR}/sample1_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPUT_DIR}/sample1_ \
     --outSAMtype BAM SortedByCoordinate

echo "¡Alineamiento completado!"
```

### Mejoras respecto a la Etapa 1

- El script es legible y tiene estructura.
- Variables parametrizadas en el encabezado.
- Se puede ejecutar (`bash rnaseq_alignment.sh`).
- Puedes versionarlo con Git.

### Problemas

- Solo procesa **una muestra** a la mano.
- Si quieres procesar otra muestra, tienes que editar el script o copiarlo.
- No hay manejo de errores.
- Si falla a mitad, no sabes dónde reempezar.

---

## Etapa 3 — Múltiples scripts, lanzamiento manual

### Descripción

Tu pipeline RNA-seq tiene ahora varios pasos: QC → Trimming → Alignment → Quantification → DEA. Cada paso es un script separado. Los lanzas **manualmente, uno después de otro**.

```
scripts/
├── 01_fastqc.sh
├── 02_trimming.sh
├── 03_alignment.sh
├── 04_featurecounts.sh
└── 05_deseq2.R
```

```bash
# Tu rutina diaria:
bash scripts/01_fastqc.sh
# [esperas 20 minutos]
bash scripts/02_trimming.sh
# [esperas 1 hora]
bash scripts/03_alignment.sh
# [esperas 3 horas, te vas a comer]
# [vuelves, ves que terminó, lanzas el siguiente]
bash scripts/04_featurecounts.sh
...
```

### Mejoras respecto a la Etapa 2

- El código está organizado en módulos lógicos.
- Cada script hace una cosa bien.
- Fácil de debuggear: sabes en qué paso falló.

### Problemas

- **Requiere supervisión constante**: alguien tiene que estar pendiente para lanzar cada paso.
- Si te vas a casa y el paso 2 falla, el paso 3 nunca empieza.
- Con 30 muestras, procesas de una en una (o duplicas y modificas scripts manualmente).
- No hay paralelización.

---

## Etapa 4 — El `.sh` maestro con envío a colas (HPC)

### Descripción

Descubres el gestor de colas de tu HPC (SLURM, SGE, LSF). Creas un script de orquestación que envía cada paso como un job, con dependencias.

```bash
#!/bin/bash
# pipeline_master.sh — Orquestador de RNA-seq para clúster

# Paso 1: FastQC para todas las muestras
JOB1=$(sbatch --array=1-30 scripts/01_fastqc.sh | awk '{print $NF}')
echo "FastQC lanzado: $JOB1"

# Paso 2: Trimming — esperar a que termine el Paso 1
JOB2=$(sbatch --dependency=afterok:$JOB1 --array=1-30 scripts/02_trimming.sh | awk '{print $NF}')
echo "Trimming lanzado: $JOB2"

# Paso 3: Alineamiento — esperar a que termine el Paso 2
JOB3=$(sbatch --dependency=afterok:$JOB2 --array=1-30 scripts/03_alignment.sh | awk '{print $NF}')
echo "Alineamiento lanzado: $JOB3"

echo "Pipeline lanzado. Jobs: $JOB1 → $JOB2 → $JOB3"
```

### Mejoras respecto a la Etapa 3

- Puedes **procesar múltiples muestras en paralelo** con job arrays.
- El pipeline corre autónomamente una vez lanzado.
- Aprovechas los recursos del clúster.
- Los pasos tienen dependencias entre sí.

### Problemas

- Frágil: si falla un job del array, los jobs dependientes no se lanzan.
- No hay logs centralizados: los logs de SLURM están dispersos por todo el sistema.
- Difícil saber qué pasó con exactamente qué muestra.
- Los paths siguen siendo específicos de tu clúster.
- Si quieres re-ejecutar solo las muestras que fallaron, tienes que hacerlo manualmente.

---

## Etapa 5 — Orquestador con logs centralizados

### Descripción

Añades un sistema de logging organizado: cada job escribe su log en una carpeta específica, con timestamps y niveles de error.

```bash
#!/bin/bash
# pipeline_master_v2.sh

LOGS_DIR="logs/$(date +%Y%m%d_%H%M%S)"
mkdir -p $LOGS_DIR

# Función para lanzar con logging
submit_with_log() {
    local STEP=$1
    local SCRIPT=$2
    local DEPENDENCY=$3
    
    local LOG_FILE="${LOGS_DIR}/${STEP}_%A_%a.log"
    
    if [ -z "$DEPENDENCY" ]; then
        JOB_ID=$(sbatch --array=1-30 \
                        --output=$LOG_FILE \
                        --error=$LOG_FILE \
                        --job-name="rnaseq_${STEP}" \
                        $SCRIPT | awk '{print $NF}')
    else
        JOB_ID=$(sbatch --array=1-30 \
                        --dependency=afterok:$DEPENDENCY \
                        --output=$LOG_FILE \
                        --error=$LOG_FILE \
                        --job-name="rnaseq_${STEP}" \
                        $SCRIPT | awk '{print $NF}')
    fi
    echo $JOB_ID
}

JOB1=$(submit_with_log "fastqc" "scripts/01_fastqc.sh")
JOB2=$(submit_with_log "trimming" "scripts/02_trimming.sh" $JOB1)
JOB3=$(submit_with_log "alignment" "scripts/03_alignment.sh" $JOB2)

echo "Pipeline lanzado. Logs en: $LOGS_DIR"
echo "Monitoriza con: squeue -u $USER"
```

```
logs/
└── 20240315_143022/
    ├── fastqc_12345_1.log
    ├── fastqc_12345_2.log
    ├── trimming_12346_1.log
    ├── alignment_12347_1.log
    └── ...
```

### Mejoras respecto a la Etapa 4

- Logs centralizados y organizados por fecha de ejecución.
- Fácil auditar qué pasó en cada paso de cada muestra.
- Nombre de job descriptivo visible en `squeue`.

### Problemas

- Sigue siendo bash: frágil, difícil de testear, difícil de mantener.
- No hay gestión de software/entornos (¿qué versión de STAR se usa?).
- No hay reproducibilidad entre sistemas.
- Re-ejecutar solo los pasos fallidos sigue requiriendo intervención manual.
- Escalar a una nueva tecnología (scRNA-seq, ATAC-seq) requiere reescribir todo.

---

## Etapa 6 — Pipeline parametrizado con metadata sheet

### Descripción

Das el salto hacia una arquitectura más robusta: **una única fuente de verdad** (el archivo de metadata) y scripts que se alimentan de ella automáticamente. El código está separado de los datos.

```
project/
├── config/
│   ├── params.config        # Parámetros globales
│   └── samples.tsv          # Metadata de todas las muestras
├── scripts/
│   ├── 01_fastqc.sh
│   ├── 02_trimming.sh
│   ├── 03_alignment.sh
│   └── 04_quantification.sh
├── run_pipeline.sh          # Punto de entrada único
└── results/
```

```tsv
# samples.tsv
sample_id   condition   fastq_R1                          fastq_R2
ctrl_1      control     /data/fastq/ctrl_1_R1.fastq.gz    /data/fastq/ctrl_1_R2.fastq.gz
ctrl_2      control     /data/fastq/ctrl_2_R1.fastq.gz    /data/fastq/ctrl_2_R2.fastq.gz
treat_1     treatment   /data/fastq/treat_1_R1.fastq.gz   /data/fastq/treat_1_R2.fastq.gz
```

```bash
# En cada script, se lee la metadata
SAMPLE_ID=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line+1{print $1}' config/samples.tsv)
FASTQ_R1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line+1{print $3}' config/samples.tsv)
```

### Mejoras respecto a la Etapa 5

- **Datos y código están separados**: añadir una muestra = añadir una fila al TSV.
- Un único punto de entrada.
- Más fácil de auditar y compartir.
- La configuración está centralizada.

### Problemas persistentes

- Aún depende completamente del entorno del clúster.
- No hay control de versiones de software integrado.
- Re-ejecuciones parciales siguen siendo manuales.
- El código sigue siendo bash con todas sus limitaciones.
- No es portable a otro HPC o a la nube sin modificaciones.

---

## Etapa 7 — Workflow Manager: el salto cuántico

### Descripción

Llegas a **Nextflow** (o Snakemake). En lugar de orquestar bash con bash, defines tu pipeline como un **grafo dirigido acíclico (DAG)** de procesos conectados por canales de datos. El sistema se encarga de la paralelización, las dependencias, los reintentos, el logging, los contenedores y la portabilidad.

```nextflow
// main.nf — Pipeline RNA-seq en Nextflow

// Definición de parámetros
params.samples_tsv  = "config/samples.tsv"
params.genome_fasta = "data/genome/Homo_sapiens.GRCh38.fa"
params.gtf          = "data/genome/Homo_sapiens.GRCh38.gtf"
params.outdir       = "results"

// Importar módulos
include { FASTQC    } from './modules/fastqc'
include { TRIMMING  } from './modules/trimming'
include { STAR_ALIGN } from './modules/star'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MULTIQC   } from './modules/multiqc'

// Canal de entrada: lee el TSV de samples
Channel
    .fromPath(params.samples_tsv)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sample_id, row.condition,
                        file(row.fastq_R1), file(row.fastq_R2)) }
    .set { ch_samples }

// Workflow principal
workflow {
    FASTQC(ch_samples)
    TRIMMING(ch_samples)
    STAR_ALIGN(TRIMMING.out.reads)
    FEATURECOUNTS(STAR_ALIGN.out.bam)
    MULTIQC(FASTQC.out.reports.collect())
}
```

### ¿Qué cambia con Nextflow?

| Característica | Pipeline en bash | Pipeline en Nextflow |
|---|---|---|
| **Paralelización** | Manual (job arrays) | Automática (por diseño) |
| **Re-ejecución parcial** | Manual | `nextflow run main.nf -resume` |
| **Portabilidad** | Específica del clúster | HPC, nube, local sin cambios |
| **Contenedores** | Configuración manual | Integrado (Docker, Singularity) |
| **Reproducibilidad** | Limitada | Total (contenedores + hash de inputs) |
| **Escalabilidad** | Complicada | Trivial (cambiar `executor`) |
| **Logs** | Dispersos | Centralizados + timeline HTML |
| **Errores** | Falla silenciosa | Informe detallado por tarea |
| **Sharing** | Difícil | `nextflow pull user/repo` |

---

## El DAG: pensar en flujos, no en pasos

Una de las revelaciones intelectuales de los workflow managers es que te obligan a pensar tu análisis como un **grafo de flujos de datos**, no como una secuencia de pasos.

```
                    ┌─────────┐
                    │ SAMPLES │
                    └────┬────┘
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
          FASTQC                TRIMMING
              │                     │
              │              ┌──────┘
              │              ▼
              │         STAR_ALIGN
              │              │
              │         FEATURECOUNTS
              │              │
              └──────┐ ┌─────┘
                     ▼ ▼
                   MULTIQC
```

Este grafo es lo que Nextflow construye y ejecuta de forma optimizada. Pasos que no tienen dependencias entre sí se ejecutan **en paralelo automáticamente**.

---

## Conclusión: La complejidad que merece la pena

Cada etapa añade complejidad, pero también añade valor:

| Etapa | Complejidad | Valor añadido |
|-------|-------------|---------------|
| 0 — Terminal | Mínima | Ninguno (no reproducible) |
| 1 — Cuaderno de comandos | Muy baja | Registro básico |
| 2 — Primer script | Baja | Reproducibilidad mínima |
| 3 — Múltiples scripts | Media | Modularidad |
| 4 — Orquestador HPC | Media-alta | Paralelización manual |
| 5 — Orquestador con logs | Alta | Auditoría |
| 6 — Pipeline parametrizado | Alta | Separación datos/código |
| 7 — Workflow Manager | Alta inicial, baja operacional | **Todo lo anterior + portabilidad + reproducibilidad total** |

> **Mensaje clave para los estudiantes:** La inversión de aprender Nextflow se amortiza en la primera vez que tienes que re-ejecutar un análisis, la primera vez que un colaborador tiene que reproducir tus resultados, y la primera vez que quieres escalar de 10 muestras a 500.

---

## Referencias y recursos

- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core: curated Nextflow pipelines](https://nf-co.re/)
- [Nextflow Training Portal](https://training.nextflow.io/)
- Di Tommaso, P. et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology 35, 316–319 (2017)
