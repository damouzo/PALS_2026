# Tips & Tricks para el Bioinformático Moderno

> *20 hacks prácticos para ser más rápido, más organizado y más efectivo en tu día a día*

---

Estos tips están pensados para bioinformáticos en etapa predoc o postdoc temprano. Muchos son pequeñas cosas que nadie te enseña explícitamente pero que, una vez las conoces, no puedes creer que hayas vivido sin ellas.

---

## 01 — Snippets personalizados en RStudio

**Qué es:** Atajos de teclado que expanden automáticamente plantillas de código en RStudio.

**Cómo configurarlo:** `Tools → Edit Code Snippets`

```
# Ejemplo de snippet personalizado:
snippet ggplot_base
	ggplot(${1:data}, aes(x = ${2:x}, y = ${3:y})) +
	  geom_point() +
	  theme_bw() +
	  labs(title = "${4:title}", x = "${2:x}", y = "${3:y}")
```

**Por qué es útil:** Escribes `ggplot_base` + Tab y tienes la plantilla completa con cursores en los puntos que necesitas modificar. Define snippets para tus funciones más repetidas: `deseq_setup`, `seurat_init`, `volcano_plot`, etc.

---

## 02 — Windows + V: el portapapeles expandido

**Qué es:** En Windows, `Win + V` abre el historial completo del portapapeles, mostrando todo lo que has copiado durante tu sesión (no solo lo último).

**Por qué es útil:** Cuántas veces has copiado algo, luego has copiado otra cosa y has perdido lo primero. Con `Win + V` tienes acceso a las últimas ~25 entradas del portapapeles. Puedes anclar las que usas frecuentemente (paths, comandos, IDs de muestras).

**Activarlo:** `Configuración → Sistema → Portapapeles → Historial del portapapeles → ON`

---

## 03 — tmux: sesiones persistentes en el servidor

**Qué es:** Un multiplexor de terminal que permite crear sesiones persistentes en el servidor remoto.

```bash
# Crear una sesión nueva
tmux new -s rnaseq_analysis

# Desconectarse (el análisis sigue corriendo)
Ctrl + B, D

# Volver a conectar
tmux attach -t rnaseq_analysis

# Listar sesiones activas
tmux ls
```

**Por qué es útil:** Si tu conexión SSH se cae o cierras el portátil, el análisis sigue corriendo. Puedes tener múltiples "ventanas" dentro de una sesión y dividir la pantalla.

---

## 04 — Alias en `.bashrc` para comandos que repites siempre

**Qué es:** Crear atajos para comandos largos o frecuentes.

```bash
# Añadir a ~/.bashrc o ~/.bash_profile:

# Monitorizar jobs en cluster
alias qme='squeue -u $USER --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R"'
alias qtop='watch -n 5 squeue -u $USER'

# Navegar a directorios frecuentes
alias proj='cd /scratch/$USER/current_project'
alias data='cd /data/shared/references'

# Git shortcuts
alias gs='git status'
alias ga='git add'
alias gc='git commit -m'

# R sin messages
alias R='R --quiet --no-save'
```

**Por qué es útil:** Cada segundo que ahorras tecleando es un segundo de más pensando. Define tus alias y actívalos con `source ~/.bashrc`.

---

## 05 — `history | grep`: encontrar comandos que ya usaste

**Qué es:** Buscar en el historial de la terminal comandos que ya ejecutaste.

```bash
# Buscar comandos que contienen "STAR"
history | grep -i star

# Buscar y re-ejecutar directamente
history | grep "nextflow run" | tail -5

# Buscar con patrón más específico
history | grep "sbatch.*rnaseq"
```

**Por qué es útil:** Cuántas veces has tenido que buscar en tus notas ese comando exacto de STAR que usaste hace 2 meses. Está en tu historial. Configura un historial grande: añade `HISTSIZE=100000` y `HISTFILESIZE=100000` a tu `.bashrc`.

---

## 06 — SSH config: aliases para tus servidores

**Qué es:** Un archivo de configuración en `~/.ssh/config` que da nombres amigables a tus conexiones SSH.

```
# ~/.ssh/config
Host hpc-queen-mary
    HostName login.hpc.qmul.ac.uk
    User qp241615
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60

Host hpc-embl
    HostName cluster.embl.de
    User tu_usuario
    ProxyJump hpc-queen-mary
```

**Por qué es útil:** En lugar de `ssh qp241615@login.hpc.qmul.ac.uk` escribes `ssh hpc-queen-mary`. Y si configuras claves SSH, sin contraseña.

---

## 07 — MultiQC: nunca veas un FastQC report por separado

**Qué es:** Herramienta que agrega los reportes de QC de múltiples herramientas (FastQC, STAR, HISAT2, featureCounts, etc.) en un único reporte HTML interactivo.

```bash
# Ejecutar MultiQC en el directorio de resultados
multiqc results/ -o results/multiqc/

# Con nombre personalizado
multiqc results/ -o results/multiqc/ --title "RNA-seq QC - Proyecto X"
```

**Por qué es útil:** Con 30 muestras, revisar 30 reportes de FastQC individualmente es imposible. MultiQC te muestra todo en un dashboard. Es estándar en cualquier pipeline de producción (nf-core lo incluye siempre). Aprende a interpretar las métricas: % de duplicados, distribución de calidad, % de alineamiento, etc.

---

## 08 — conda/mamba environments: adiós a "pero en mi máquina funciona"

**Qué es:** Entornos virtuales que encapsulan todas las dependencias de software de un análisis.

```bash
# Crear entorno para RNA-seq
conda create -n rnaseq python=3.10 star=2.7.11 samtools=1.18 fastqc

# Activar
conda activate rnaseq

# Exportar para reproducibilidad
conda env export > environment.yml

# Recrear en otro sistema
conda env create -f environment.yml
```

**Por qué es útil:** Cada análisis en su propio entorno. No más conflictos de versiones entre proyectos. El archivo `environment.yml` es parte de la reproducibilidad del análisis: inclúyelo en tu repositorio.

**Tip adicional:** usa `mamba` en lugar de `conda` para instalaciones más rápidas.

---

## 09 — Quarto / RMarkdown: el cuaderno de laboratorio computacional

**Qué es:** Sistema de documentos que mezcla texto narrativo, código R o Python, y resultados (gráficos, tablas) en un único documento HTML, PDF o Word.

```r
---
title: "Análisis DEG - Proyecto KO"
author: "Tu Nombre"
date: "`r Sys.Date()`"
output: html_document
params:
  condition: "KO_vs_WT"
  fdr_threshold: 0.05
---

## Resultados DEseq2

```{r}
# El código y los resultados van juntos
res <- results(dds, contrast = c("condition", params$condition))
plotMA(res, alpha = params$fdr_threshold)
```
```

**Por qué es útil:** Tu análisis es su propia documentación. Puedes reproducirlo exactamente. Fácil de compartir con tu supervisor. Los parámetros parametrizados permiten re-ejecutar el reporte con distintas condiciones sin tocar el código.

---

## 10 — Git: versionar el código, no los datos

**Qué es:** Control de versiones para tu código. El historial completo de cambios en tu pipeline, scripts y análisis.

```bash
# Flujo de trabajo básico
git init
git add scripts/ analysis/ nextflow.config
git commit -m "Add DESeq2 analysis with batch correction"

# Ver qué ha cambiado
git diff

# Volver a una versión anterior del análisis
git log --oneline
git checkout abc1234 -- scripts/deseq2_analysis.R
```

**Por qué es útil:** "¿Qué cambié exactamente en el script para obtener estos resultados?" — Git lo sabe. "Quiero probar este cambio sin romper lo que ya funciona" — las ramas de Git están hechas para eso.

**Regla de oro:** Los datos van en almacenamiento de datos (no en Git), el código va en Git. El archivo `.gitignore` es tu amigo.

---

## 11 — `ln -s`: enlaces simbólicos para no duplicar datos

**Qué es:** Un acceso directo del sistema de archivos que apunta a otro archivo o directorio.

```bash
# Crear un enlace simbólico a los datos originales
ln -s /data/shared/references/GRCh38/genome.fa data/genome.fa

# En lugar de copiar 30 GB, el enlace ocupa bytes
ls -la data/genome.fa
# lrwxrwxrwx 1 user group 45 Jan 1 12:00 data/genome.fa -> /data/shared/references/...
```

**Por qué es útil:** No duplicas datos pesados (genomas, índices). Tu estructura de proyecto queda limpia. Si los datos originales se actualizan, el enlace apunta a la versión correcta.

---

## 12 — `rsync`: copiar datos de forma inteligente

**Qué es:** Herramienta de copia que solo transfiere lo que ha cambiado, con soporte de reanudación.

```bash
# Copiar resultados al servidor remoto (solo lo nuevo)
rsync -avzh --progress results/ servidor:/home/user/results/

# Con compresión y checkedsum
rsync -avzc --exclude='work/' results/ backup:/archive/project_x/

# Vista previa sin transferir
rsync -avzhn results/ servidor:/home/user/results/
```

**Por qué es útil:** `scp` copia todo siempre. `rsync` es inteligente: si una transferencia de 50 GB se corta a mitad, la reanuda desde donde la dejó.

---

## 13 — `watch`: monitorizar comandos de forma continua

**Qué es:** Re-ejecuta un comando cada N segundos mostrando la salida actualizada.

```bash
# Ver el estado de tus jobs en SLURM actualizado cada 5 segundos
watch -n 5 squeue -u $USER

# Ver el crecimiento de un directorio de resultados
watch -n 10 'ls -la results/ | wc -l'

# Ver si un proceso de Nextflow está progresando
watch -n 30 'cat .nextflow.log | tail -20'
```

**Por qué es útil:** En lugar de ejecutar `squeue` manualmente cada 2 minutos, `watch` lo hace por ti y resalta qué ha cambiado.

---

## 14 — `.Rprofile`: personalizar tu entorno de R

**Qué es:** Un script de R que se ejecuta automáticamente al iniciar R o RStudio.

```r
# ~/.Rprofile
options(
  stringsAsFactors = FALSE,  # Jamás strings como factors automáticamente
  warn = 1,                  # Warnings inmediatos (no al final)
  scipen = 999,              # Sin notación científica
  repos = c(CRAN = "https://cloud.r-project.org")
)

# Función para cargar tus paquetes frecuentes rápido
.my_libs <- function() {
  library(tidyverse)
  library(ggplot2)
  library(DESeq2)
  cat("Paquetes cargados!\n")
}

# Autocompletar con Tab en la consola
utils::rc.settings(ipck = TRUE)

cat("R listo!\n")
```

**Por qué es útil:** Cada vez que abres R tienes tu entorno configurado a tu gusto. La función `.my_libs()` carga tus paquetes más usados con un solo comando.

---

## 15 — `md5sum`: verificar integridad de datos

**Qué es:** Genera o comprueba un hash criptográfico de un archivo para verificar que no se ha corrompido.

```bash
# Generar checksums de los FASTQs recibidos
md5sum data/fastq/*.fastq.gz > data/fastq/md5sums.txt

# Verificar integridad después de una transferencia
md5sum -c data/fastq/md5sums.txt
# sample1_R1.fastq.gz: OK
# sample1_R2.fastq.gz: OK

# Comparar con los checksums del proveedor de secuenciación
diff md5sums_recibidos.txt md5sums.txt
```

**Por qué es útil:** Un FASTQ corrompido puede dar resultados de alineamiento inesperados sin dar error obvio. Verifica siempre la integridad de los datos, especialmente después de transferencias largas.

---

## 16 — `htop` / `top`: ver qué está consumiendo tu servidor

**Qué es:** Monitores interactivos de los procesos y recursos del sistema.

```bash
# htop (más visual, si está instalado)
htop

# Filtrar por usuario
htop -u $USER

# top clásico
top
```

**Por qué es útil:** Antes de lanzar un análisis intensivo, comprueba que el servidor no está saturado. Después de un análisis, puedes ver exactamente cuánta RAM y CPU consumió para dimensionar correctamente tus jobs de SLURM.

---

## 17 — Variables de entorno en `.bashrc` para paths frecuentes

**Qué es:** Variables que apuntan a rutas frecuentes, disponibles en toda la terminal.

```bash
# ~/.bashrc
export REF_GENOME="/data/references/GRCh38"
export MY_PROJ="/scratch/$USER/current_project"
export SCRIPTS_DIR="$MY_PROJ/scripts"

# Ahora puedes usarlas en cualquier comando
STAR --genomeDir $REF_GENOME/star_index ...
cd $MY_PROJ
bash $SCRIPTS_DIR/01_fastqc.sh
```

**Por qué es útil:** Cuando cambias de proyecto o el path del genoma de referencia cambia, solo tienes que actualizar un lugar. Todos los scripts que usan la variable se actualizan automáticamente.

---

## 18 — VSCode Remote SSH: tu IDE en el servidor remoto

**Qué es:** La extensión Remote - SSH de VS Code te permite editar archivos directamente en el servidor remoto como si estuvieras en local.

**Cómo configurarlo:**
1. Instalar extensión "Remote - SSH" en VS Code
2. `Ctrl + Shift + P` → "Remote-SSH: Connect to Host..."
3. Introducir `usuario@servidor.hpc.ac.uk`
4. VS Code se conecta y puedes editar, navegar y usar el terminal integrado en el servidor

**Por qué es útil:** Adiós a editar scripts con `nano` o `vim` en el servidor, y a la pesadilla de transferir archivos cada vez que haces un cambio. Con la extensión de Nextflow para VS Code, además tienes syntax highlighting y autocompletado.

---

## 19 — `grep` con contexto: buscar en logs de forma inteligente

**Qué es:** Opciones de `grep` para ver contexto alrededor de la línea que buscas.

```bash
# Ver 3 líneas antes y después del error
grep -B 3 -A 5 "ERROR" .nextflow.log

# Buscar recursivamente en todos los logs del directorio work/
grep -r "OutOfMemoryError" work/

# Contar cuántas muestras fallaron
grep -c "ERROR" work/*/command.log

# Ver el error de una task específica
cat work/82/4574828abc123/.command.err
```

**Por qué es útil:** Los errores en bioinformática raramente están solos. Ver las líneas de contexto ayuda a entender qué pasó exactamente. Esencial para debuggear pipelines de Nextflow.

---

## 20 — GitHub Codespaces / Gitpod: tu entorno de desarrollo en la nube

**Qué es:** Entornos de desarrollo completamente configurados que corren en la nube y se abren en el navegador.

```
# Para usar el Codespace oficial de Nextflow training:
1. Ir a https://training.nextflow.io
2. Hacer clic en "Open in GitHub Codespaces"
3. Tienes un entorno con Nextflow, Docker y todo configurado en 2 minutos
```

**Por qué es útil:** No necesitas instalar nada. Ideal para talleres, cursos o para probar algo rápidamente sin tocar tu entorno local. El Codespace oficial de nf-core/training es especialmente útil para aprender Nextflow.

---

## Bonus: Los 5 comandos que más tiempo me han ahorrado

```bash
# 1. Re-ejecutar el último comando como root (o con sudo)
sudo !!

# 2. Volver al directorio anterior
cd -

# 3. Repetir el último argumento del comando anterior
cp file_largo.txt !!:$   # equivale a cp file_largo.txt (último arg)

# 4. Expandir el último argumento con Alt+.
# Escribe Alt+. y aparece el último argumento del comando previo

# 5. Ejecutar un comando en segundo plano sin que muera al cerrar la terminal
nohup bash mi_script.sh > mi_script.log 2>&1 &
```

---

## Recursos para seguir descubriendo

- [The Missing Semester of Your CS Education (MIT)](https://missing.csail.mit.edu/) — El curso que te enseña todo lo que los cursos de programación no enseñan
- [Bioinformatics one-liners](https://github.com/stephenturner/oneliners) — Colección de one-liners para manipulación de datos
- [Cheat.sh](https://cheat.sh/) — Documentación rápida de cualquier comando desde la terminal (`curl cheat.sh/rsync`)
- [explainshell.com](https://explainshell.com/) — Pega cualquier comando de bash y te explica cada parte
