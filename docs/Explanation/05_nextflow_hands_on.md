# Nextflow Hands-On: Tus Primeros Pasos

> *Guía práctica para la sesión hands-on del seminario*  
> **Basado en Nextflow DSL2 — versión 25.x**  
> **Referencia oficial:** [nextflow.io/docs/latest/your-first-script.html](https://nextflow.io/docs/latest/your-first-script.html)

---

## Antes de empezar

### Requisitos

```bash
# Verificar Java (Nextflow requiere Java >= 11)
java -version

# Verificar que Nextflow está instalado
nextflow -version
# Debería mostrar algo como: nextflow version 25.10.0...
```

### ¿No tienes Nextflow instalado?

```bash
# Instalación en una línea (Linux/macOS)
curl -s https://get.nextflow.io | bash

# Mover a un directorio en tu PATH
mv nextflow ~/bin/
```

### Alternativa: GitHub Codespaces (sin instalar nada)

Si tienes cuenta de GitHub, puedes usar el entorno oficial de training:
1. Ve a [github.com/nextflow-io/training](https://github.com/nextflow-io/training)
2. Haz clic en `Code → Codespaces → Create codespace`
3. Tendrás un entorno completamente configurado en el navegador

---

## Parte 1 — Hello World: tu primer pipeline

### 1.1 Crear el archivo del pipeline

Crea un directorio de trabajo y el archivo principal:

```bash
mkdir nextflow_hands_on && cd nextflow_hands_on
```

Crea el archivo `main.nf` con el siguiente contenido:

```nextflow
// main.nf
// ============================================================
// Pipeline Hello World: split y conversión a mayúsculas
// ============================================================

// Parámetro de entrada con valor por defecto
params.str = "Hello world!"

// ── Proceso 1: split ──────────────────────────────────────
// Divide una cadena de texto en chunks de 6 caracteres
process split {
    input:
    val x        // Recibe un valor (string)

    output:
    path 'chunk_*'  // Produce archivos con prefijo 'chunk_'

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// ── Proceso 2: convert_to_upper ───────────────────────────
// Convierte el contenido de un archivo a mayúsculas
process convert_to_upper {
    tag "$y"  // Etiqueta para identificar la task en los logs

    input:
    path y    // Recibe un archivo

    output:
    path 'upper_*'  // Produce archivos con prefijo 'upper_'

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}

// ── Bloque de output ──────────────────────────────────────
// Define cómo y dónde se publican los resultados
output {
    lower {
        path 'lower'
    }
    upper {
        path 'upper'
    }
}

// ── Workflow principal ────────────────────────────────────
workflow {
    main:
    // 1. Crear un canal con el parámetro de entrada
    ch_str = channel.of(params.str)

    // 2. Ejecutar split y obtener los chunks
    ch_chunks = split(ch_str)

    // 3. Convertir a mayúsculas (flatten: un archivo por tarea)
    ch_upper = convert_to_upper(ch_chunks.flatten())

    publish:
    lower = ch_chunks.flatten()  // Publicar chunks originales
    upper = ch_upper             // Publicar resultados en mayúsculas
}
```

### 1.2 Ejecutar el pipeline

```bash
nextflow run main.nf
```

**Salida esperada:**

```
 N E X T F L O W   ~  version 25.10.0

Launching `main.nf` [big_wegener] DSL2 - revision: 13a41a8946

executor >  local (3)
[82/457482] split (1)                   | 1 of 1 ✔
[2f/056a98] convert_to_upper (chunk_aa) | 2 of 2 ✔
```

**Explorar los resultados:**

```bash
# Ver la estructura de resultados
ls -R results/

# Ver el contenido
cat results/lower/chunk_aa
cat results/upper/upper_chunk_aa
```

### 1.3 Entender el directorio `work/`

Nextflow aísla cada task en su propio subdirectorio dentro de `work/`:

```bash
# Ver la estructura del directorio work
ls work/

# Explorar una task específica (el hash viene del log)
ls work/82/457482.../
```

Dentro de cada directorio de task encontrarás:
- `.command.sh` — el script exacto que se ejecutó
- `.command.log` — el output completo (stdout + stderr)
- `.command.err` — stderr separado
- `.exitcode` — el código de salida (0 = éxito)
- Los archivos de entrada (como symlinks) y de salida

```bash
# Ver el script que ejecutó Nextflow
cat work/82/457482*/.command.sh

# Ver el log de la task
cat work/82/457482*/.command.log
```

> **Concepto clave:** Cada task se ejecuta en su propio directorio aislado. Esto es lo que hace que Nextflow pueda reconstruir exactamente qué se ejecutó, con qué inputs.

---

## Parte 2 — Resume: la magia del caché

### 2.1 Modificar el pipeline y re-ejecutar

Modifica el proceso `convert_to_upper` para que **invierta** el texto en lugar de convertirlo a mayúsculas:

```nextflow
// Reemplaza el proceso convert_to_upper por este:
process convert_to_upper {
    tag "$y"

    input:
    path y

    output:
    path 'upper_*'

    script:
    """
    rev $y > upper_${y}
    """
}
```

Ahora ejecuta con `-resume`:

```bash
nextflow run main.nf -resume
```

**Salida esperada:**

```
N E X T F L O W   ~  version 25.10.0

Launching `main.nf` [furious_curie] DSL2 - revision: 5490f13c43

executor >  local (2)
[82/457482] split (1)                   | 1 of 1, cached: 1 ✔
[02/9db40b] convert_to_upper (chunk_aa) | 2 of 2 ✔
```

**Observaciones:**
- El proceso `split` muestra `cached: 1` — **no se re-ejecutó**
- El proceso `convert_to_upper` sí se re-ejecutó porque el código cambió
- El output del proceso `split` se reutilizó del caché

> **Concepto clave:** Nextflow calcula un hash de (código del proceso + inputs + configuración). Si el hash no cambia, la task se cachea. Esto es fundamentalmente diferente al sistema de timestamps de Snakemake o Make.

### 2.2 Cambiar el parámetro de entrada

```bash
nextflow run main.nf --str 'Bonjour le monde'
```

Observa cómo ahora hay más chunks porque el string es más largo:

```
executor >  local (4)
[55/a3a700] split (1)                   | 1 of 1 ✔
[f4/af5ddd] convert_to_upper (chunk_aa) | 3 of 3 ✔
```

---

## Parte 3 — De Hello World a Bioinformática

### 3.1 Un pipeline RNA-seq mínimo

Vamos a estructurar un pipeline que procese samples de RNA-seq. Crea un nuevo directorio:

```bash
mkdir rnaseq_mini && cd rnaseq_mini
```

Crea la estructura del proyecto:

```bash
mkdir -p data/fastq modules conf
```

Crea un archivo de metadata de samples (`conf/samples.tsv`):

```tsv
sample_id	condition	fastq_r1	fastq_r2
ctrl_1	control	data/fastq/ctrl_1_R1.fastq.gz	data/fastq/ctrl_1_R2.fastq.gz
ctrl_2	control	data/fastq/ctrl_2_R1.fastq.gz	data/fastq/ctrl_2_R2.fastq.gz
treat_1	treatment	data/fastq/treat_1_R1.fastq.gz	data/fastq/treat_1_R2.fastq.gz
treat_2	treatment	data/fastq/treat_2_R1.fastq.gz	data/fastq/treat_2_R2.fastq.gz
```

### 3.2 El pipeline principal (`main.nf`)

```nextflow
// main.nf — Pipeline RNA-seq mínimo en Nextflow DSL2
nextflow.enable.dsl=2

// ── Parámetros ────────────────────────────────────────────
params.samples_tsv  = "conf/samples.tsv"
params.outdir       = "results"
params.genome_index = "/path/to/genome/index"  // Cambiar a tu path

// ── Importar módulos ──────────────────────────────────────
include { FASTQC         } from './modules/fastqc'
include { TRIMGALORE     } from './modules/trimgalore'
include { STAR_GENOMEGENERATE } from './modules/star_index'
include { STAR_ALIGN     } from './modules/star_align'

// ── Workflow principal ────────────────────────────────────
workflow {

    // 1. Crear canal de samples desde el TSV
    ch_samples = Channel
        .fromPath(params.samples_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            // Devuelve un tuple: [meta, [fastq_R1, fastq_R2]]
            def meta = [id: row.sample_id, condition: row.condition]
            def reads = [file(row.fastq_r1), file(row.fastq_r2)]
            return [meta, reads]
        }

    // 2. QC de las lecturas crudas
    FASTQC(ch_samples)

    // 3. Trimming de adaptadores
    TRIMGALORE(ch_samples)

    // 4. Alineamiento con STAR
    STAR_ALIGN(
        TRIMGALORE.out.reads,  // Output del proceso anterior
        params.genome_index
    )
}
```

### 3.3 Un módulo simple: FASTQC (`modules/fastqc.nf`)

```nextflow
// modules/fastqc.nf
// Módulo para ejecutar FastQC sobre reads paired-end

process FASTQC {
    tag "${meta.id}"  // El tag identifica la task en el log

    // Dónde se publican los resultados
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    // Software (conda o container)
    conda "bioconda::fastqc=0.12.1"
    // container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)  // meta map + lista de archivos

    output:
    tuple val(meta), path("*.{html,zip}"), emit: reports

    script:
    """
    fastqc \\
        --threads 4 \\
        --outdir . \\
        ${reads}
    """
}
```

> **Concepto: el meta map**
> El patrón `tuple val(meta), path(reads)` es el idioma estándar de nf-core.
> `meta` es un mapa (diccionario) con información del sample: `[id: "ctrl_1", condition: "control"]`.
> Esto permite propagar metadatos a lo largo de todo el pipeline sin perder el contexto de cada sample.

### 3.4 Configuración (`nextflow.config`)

```groovy
// nextflow.config — Configuración del pipeline

// Parámetros por defecto
params {
    outdir          = "results"
    samples_tsv     = "conf/samples.tsv"
    genome_index    = null
    max_cpus        = 16
    max_memory      = '64.GB'
    max_time        = '24.h'
}

// Perfiles de ejecución
profiles {
    local {
        executor.name = 'local'
        executor.cpus = 8
        executor.memory = '32 GB'
    }

    slurm {
        executor.name = 'slurm'
        executor.queueSize = 50
        process.queue = 'long'
        process.clusterOptions = '--account=your_account'
    }

    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
    }

    test {
        params.samples_tsv  = "test/samples_test.tsv"
        params.genome_index = "test/genome_index"
        executor.name       = 'local'
    }
}

// Configuración por defecto de procesos
process {
    // Recursos por defecto
    cpus   = 2
    memory = '8.GB'
    time   = '4.h'

    // Labels para asignar recursos específicos
    withLabel: 'process_low' {
        cpus   = 2
        memory = '8.GB'
        time   = '4.h'
    }
    withLabel: 'process_medium' {
        cpus   = 8
        memory = '32.GB'
        time   = '8.h'
    }
    withLabel: 'process_high' {
        cpus   = 16
        memory = '64.GB'
        time   = '24.h'
    }

    // Reintentar con más memoria si falla
    errorStrategy = { task.exitStatus in [137,140] ? 'retry' : 'finish' }
    maxRetries    = 2
    memory        = { task.attempt <= 2 ? task.memory * task.attempt : task.memory }
}

// Reportes automáticos
timeline {
    enabled = true
    file    = "${params.outdir}/pipeline_info/timeline.html"
}
report {
    enabled = true
    file    = "${params.outdir}/pipeline_info/report.html"
}
dag {
    enabled = true
    file    = "${params.outdir}/pipeline_info/dag.html"
}
```

### 3.5 Cómo ejecutar el pipeline

```bash
# En local (con conda)
nextflow run main.nf -profile local

# En HPC con SLURM y Singularity
nextflow run main.nf -profile slurm,singularity

# Re-ejecutar desde donde lo dejaste
nextflow run main.nf -resume -profile slurm,singularity

# Con parámetros personalizados
nextflow run main.nf \
    --samples_tsv conf/my_samples.tsv \
    --outdir /scratch/user/results \
    -profile slurm,singularity
```

---

## Parte 4 — Usar un pipeline nf-core

### 4.1 Ejecutar nf-core/rnaseq en modo test

El perfil `test` usa un dataset pequeño incluido en el pipeline para validar que funciona:

```bash
# Asegurarse de que nextflow y Docker/Singularity están disponibles
nextflow -version
docker --version  # o singularity --version

# Ejecutar el test de nf-core/rnaseq
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir results_test
```

Esto descarga el pipeline directamente desde GitHub y lo ejecuta con datos de test. El proceso toma ~10-15 minutos.

### 4.2 Preparar un análisis real con nf-core/rnaseq

**El samplesheet** (formato CSV, requerido por nf-core/rnaseq):

```csv
sample,fastq_1,fastq_2,strandedness
ctrl_rep1,/data/fastq/ctrl_1_R1.fastq.gz,/data/fastq/ctrl_1_R2.fastq.gz,auto
ctrl_rep2,/data/fastq/ctrl_2_R1.fastq.gz,/data/fastq/ctrl_2_R2.fastq.gz,auto
treat_rep1,/data/fastq/treat_1_R1.fastq.gz,/data/fastq/treat_1_R2.fastq.gz,auto
treat_rep2,/data/fastq/treat_2_R1.fastq.gz,/data/fastq/treat_2_R2.fastq.gz,auto
```

**Ejecutar el análisis completo:**

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile singularity \
    -resume
```

### 4.3 Explorar el output de nf-core/rnaseq

```
results/
├── fastqc/          # QC de reads crudos
├── trimgalore/      # Reads después de trimming
├── star_salmon/     # Alineamientos BAM y cuantificación
├── salmon/          # Cuantificación con Salmon
├── deseq2_qc/       # QC de expresión con DESeq2
├── multiqc/         # Reporte QC integrado (¡empieza aquí!)
└── pipeline_info/   # Timeline, DAG, reporte de ejecución
```

---

## Parte 5 — Tips de debugging

### 5.1 Ver los logs de una task fallida

```bash
# El log de Nextflow
cat .nextflow.log | grep "ERROR" -A 10

# Encontrar la task que falló
nextflow log | head  # Ver las últimas ejecuciones

# Ver los logs de error de tasks específicas
cat work/ab/12345*/.command.err

# El script exacto que se ejecutó
cat work/ab/12345*/.command.sh
```

### 5.2 Entrar en el directorio de una task para debuggear

```bash
# Ir al directorio de la task
cd work/ab/12345abc456.../

# Ver qué hay
ls -la

# Re-ejecutar el script manualmente para ver el error
bash .command.sh
```

### 5.3 Comandos de diagnóstico de Nextflow

```bash
# Ver el historial de ejecuciones
nextflow log

# Ver los detalles de la última ejecución
nextflow log last

# Ver todas las tasks de una ejecución específica
nextflow log last -f name,status,exit,work_dir

# Ver solo las tasks que fallaron
nextflow log last -f name,status,exit,work_dir -filter 'status == "FAILED"'
```

### 5.4 Errores comunes y soluciones

| Error | Causa probable | Solución |
|-------|----------------|----------|
| `No such file or directory` en el script | Path incorrecto en el input | Verificar los paths en el canal de entrada |
| `command not found` | Software no instalado en el entorno | Verificar el entorno conda o la imagen Docker |
| `Exit code 137` | Proceso killed por Out Of Memory | Aumentar la memoria del proceso con un label |
| `WARN: Not a valid cache. Will retry from scratch` | Cache corrupted | Eliminar `.nextflow/` y re-ejecutar |
| `Channel value was already consumed` | Re-uso de canal sin `view()` | Usar `.tap {}` o reestructurar el workflow |

---

## Resumen de comandos clave

```bash
# Ejecutar
nextflow run main.nf

# Re-ejecutar desde caché
nextflow run main.nf -resume

# Ejecutar con parámetro personalizado
nextflow run main.nf --str "Mi texto"

# Ejecutar con perfil
nextflow run main.nf -profile slurm,singularity

# Ver historial
nextflow log

# Ver tareas de la última ejecución
nextflow log last -f name,status,exit

# Limpiar el directorio work (¡cuidado!)
nextflow clean -f

# Ver versión
nextflow -version

# Actualizar Nextflow
nextflow self-update
```

---

## Estructura recomendada de un proyecto Nextflow

```
mi_pipeline/
├── main.nf                  # Workflow principal
├── nextflow.config          # Configuración y perfiles
├── nextflow_schema.json     # Schema de parámetros (nf-core standard)
├── README.md                # Documentación del pipeline
├── CHANGELOG.md             # Historial de cambios
│
├── conf/                    # Configuraciones específicas
│   ├── base.config          # Recursos por defecto
│   ├── igenomes.config      # Genomas de referencia
│   └── test.config          # Configuración para tests
│
├── modules/                 # Módulos de procesos
│   ├── local/               # Módulos propios
│   │   └── fastqc/
│   │       └── main.nf
│   └── nf-core/             # Módulos instalados de nf-core
│       └── fastqc/
│           └── main.nf
│
├── subworkflows/            # Sub-workflows reutilizables
│   └── qc_reads/
│       └── main.nf
│
├── bin/                     # Scripts auxiliares
│   └── parse_samples.py
│
└── test/                    # Datos y configuración de test
    ├── samples_test.tsv
    └── data/
```

---

## Próximos pasos

Una vez que te sientas cómodo con lo básico:

1. **[Hello Nextflow](https://training.nextflow.io/latest/hello_nextflow/)** — Tutorial oficial completo (Channels, Modules, Containers, Config)
2. **[Nextflow for RNAseq](https://training.nextflow.io/latest/nf4_science/rnaseq/)** — Aplicación directa a tu análisis
3. **[nf-core modules](https://nf-co.re/modules)** — Biblioteca de módulos listos para usar
4. **[nf-core/tools](https://nf-co.re/tools)** — CLI para crear pipelines con la estructura estándar nf-core
5. **[Side Quests: Debugging Workflows](https://training.nextflow.io/latest/side_quests/debugging/)** — Cómo debuggear como un experto

---

## Recursos de referencia

| Recurso | URL |
|---------|-----|
| Documentación oficial | [nextflow.io/docs/latest/](https://nextflow.io/docs/latest/) |
| Portal de training | [training.nextflow.io](https://training.nextflow.io/) |
| nf-core pipelines | [nf-co.re/pipelines](https://nf-co.re/pipelines) |
| nf-core modules | [nf-co.re/modules](https://nf-co.re/modules) |
| Foro de la comunidad | [community.seqera.io](https://community.seqera.io/) |
| nf-core Slack | [nf-co.re/join/slack](https://nf-co.re/join/slack) |
| Canal de YouTube | [youtube.com/@Nextflow](https://www.youtube.com/@Nextflow) |
