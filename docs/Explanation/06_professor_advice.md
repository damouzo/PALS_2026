# Consejo del Catedrático: Revisión del Seminario de Workflow Management

> *¿Qué diría un catedrático senior en bioinformática computacional si le presentas este plan de clase?*

---

## El escenario

Imagina que le presentas tu plan de seminario a un catedrático con 20 años de experiencia en bioinformática computacional, que ha formado a decenas de estudiantes de doctorado, ha publicado en Nature Methods, ha colaborado con el consorcio de nf-core y ha visto evolucionar el campo desde los scripts de Perl hasta los pipelines de producción modern. Esta sería su respuesta.

---

## "Me parece una propuesta sólida y bien estructurada. Aquí están mis observaciones."

---

### Lo que está muy bien planteado

**El arco narrativo de la evolución es brillante.**  
Empezar con "así es como todos empezamos" crea una conexión inmediata con la audiencia. Los predocs se reconocerán en las primeras etapas, y eso genera la motivación intrínseca para querer llegar a las últimas. Es mucho más efectivo pedagógicamente que empezar directamente con "Nextflow es un workflow manager que...". *No toques eso.*

**Centrarte en Nextflow en lugar de dar una panorámica igualitaria de todas las herramientas es la decisión correcta.**  
Un seminario de 3-4 horas en el que intentas cubrir Snakemake, Nextflow, WDL y CWL con equidad no sirve para nada. La profundidad en una herramienta siempre supera a la amplitud superficial. La audiencia saldrá con algo accionable, no con confusión.

**Los tips y tricks son un acierto.**  
En mi experiencia, es lo que mejor recuerdan los estudiantes semanas después. Cosas concretas, pequeñas victorias inmediatas. Bien como cierre o como "bonus track" que da energía.

---

### Lo que añadiría o cambiaría

#### 1. El "por qué importa" necesita más fuerza emocional

Tu argumento sobre reproducibilidad es correcto, pero en esta audiencia (predocs que están pensando en sobrevivir su próxima reunión de supervisión), el argumento más poderoso no es "la ciencia debe ser reproducible". Es mucho más cercano:

> *"Dentro de 6 meses, tu supervisor te va a pedir que repitas ese análisis con 10 muestras más. ¿Qué vas a hacer?"*

> *"Acabas de recibir los comentarios del revisor número 2, que pide que cambies el umbral de FDR y recalcules todos los resultados. ¿Cuánto tiempo te va a costar?"*

Estos escenarios son concretos, son universales y crean urgencia real. Considera abrir con uno de ellos en lugar de con la definición académica de reproducibilidad.

#### 2. Falta abordar el miedo a romper algo en producción

Uno de los mayores inhibidores para adoptar herramientas nuevas es el miedo: *"¿y si lo rompo?"*, *"¿y si mi pipeline de siempre deja de funcionar mientras aprendo esto?"*, *"¿y si mi supervisor me pregunta por qué tardé más?"*. 

Dedicaría 5 minutos a hablar explícitamente de esto. El mensaje: **Nextflow no reemplaza tus scripts de hoy para mañana. Es una habilidad que construyes en paralelo.** Empieza con el pipeline más sencillo de tu análisis, no con el más complejo.

#### 3. Habla explícitamente de cómo vender Nextflow a tu supervisor

Los predocs no solo tienen que convencerse a sí mismos: tienen que convencer a su PI. Muchos PIs tienen la percepción de que aprender Nextflow es "perder tiempo" que debería dedicarse a producir resultados. 

Sugiero añadir un slide o una sección titulada *"Cómo hablar de esto con tu supervisor"* con los argumentos concretos:
- "El mismo análisis me tardó 2 días en lanzar manualmente; con Nextflow tarda 20 minutos y puedo re-ejecutarlo en cualquier momento"
- "Cuando enviemos el paper, el revisor podrá ejecutar exactamente el mismo pipeline con los mismos resultados"
- "Cuando nuestro colaborador del Broad quiera reproducir el análisis, tenemos una línea de comando"

#### 4. La sección de contenedores merece más protagonismo

Los contenedores (Docker/Singularity) son, junto con el DAG, el cambio más profundo que trae un workflow manager. Y sin embargo, muchos bioinformáticos los evitan porque los ven como "cosa de ingenieros de software". 

Yo añadiría un capítulo específico dedicado a:
- Qué es una imagen de contenedor en términos simples (receta de cocina reproducible para tu software)
- Cómo usar las imágenes de Biocontainers (prácticamente cualquier herramienta bioinformática ya tiene su imagen)
- Cómo especificar el contenedor en un proceso de Nextflow
- Por qué `singularity` es el estándar en HPC y `docker` en local/cloud

Sin contenedores, la reproducibilidad de Nextflow es parcial.

#### 5. Considera añadir un bloque sobre el ciclo de vida de un análisis real

Hay un gap entre "entiendo Nextflow conceptualmente" y "sé cómo estructuro mi análisis de RNA-seq en un pipeline de Nextflow desde cero". Ese gap es enorme y desalienta a la gente.

Propongo mostrar un ejemplo real, aunque sea simplificado, que siga el ciclo completo:

```
Datos crudos → QC → Trimming → Alineamiento → Cuantificación → DEA → Visualización
```

No tiene que ser funcional al 100%, pero mostrar cómo ese `main.nf` se estructura, cómo los outputs de un proceso se convierten en inputs del siguiente, y dónde encajan los módulos de nf-core. Esto es lo que convierte el Hello World abstracto en algo aplicable a su trabajo.

#### 6. Habla de los errores que todos cometen al empezar

Años de experiencia enseñando bioinformática me han enseñado que los estudiantes cometen los mismos errores al iniciarse en Nextflow. Sería valioso incluir aunque sea una lista corta:

- **Hardcodear paths absolutos en el script de proceso**: un proceso que tiene `/scratch/user123/data/...` en el código no es portátil. Los inputs deben venir siempre de los canales.
- **Re-usar canales**: en Nextflow, un canal se consume. Si intentas usar el mismo canal dos veces, el segundo uso estará vacío. Solución: `.tap()` o redefinir.
- **No usar `tag`**: sin tag, en un pipeline con 100 muestras todos los logs dicen `FASTQC (1)`. Con `tag "${meta.id}"` dices `FASTQC (ctrl_sample_A)`.
- **Publicar outputs con `mode: 'move'`**: si mueves los archivos de `work/` a `results/`, pierdes el caché para `-resume`. Usa siempre `mode: 'copy'` o `mode: 'symlink'`.
- **Ignorar el directorio `work/` hasta que algo falla**: enséñales desde el principio a usarlo como herramienta de debugging, no solo como directorio temporal.

#### 7. nf-core antes que escribir desde cero

Este es quizás mi consejo más importante y el que cambia completamente el enfoque del primer día con Nextflow:

> **El primer Nextflow que deberías ejecutar es `nextflow run nf-core/rnaseq -profile test,docker`. No es `main.nf` desde cero.**

El orden que yo recomendaría:
1. **Primero, correr un pipeline nf-core completo** — entender qué hace, explorar el output, aprender qué es un samplesheet
2. **Segundo, inspeccionar el código del pipeline nf-core** — ver cómo está estructurado, qué módulos usa, cómo está configurado
3. **Tercero, Hello World** — ahora que has visto el destino final, construir desde los fundamentos tiene más sentido
4. **Cuarto, módulos propios** — adaptar módulos de nf-core a tus necesidades específicas

Este orden invierte la curva de frustración: empiezas con una victoria grande (tengo el RNA-seq completo funcionando) en lugar de con un ejercicio abstracto (convirtiendo texto a mayúsculas).

---

### Lo que me parece importante y estás omitiendo

#### Control de versiones y trazabilidad del análisis

El seminario habla de Nextflow pero no aborda suficientemente la pregunta: *"¿Cómo sé exactamente qué versión del pipeline produzco estos resultados?"*

Esto incluye:
- **Git tagging**: etiquetar el commit exacto que produjo los resultados de una figura en un paper
- **`nextflow.config` con versiones pinned**: no `nextflow run nf-core/rnaseq` sino `nextflow run nf-core/rnaseq -r 3.14.0`
- **El `pipeline_info/` que genera nf-core**: contiene el exact software version report
- **El manifest en el `nextflow.config`** para documentar versión del pipeline

Esta trazabilidad es lo que te salva cuando el reviewer pregunta "¿qué versión de STAR usaste?"

#### Gestión de datos: el problema más grande que nadie controla bien

Nextflow gestiona el análisis, pero los **datos** son una responsabilidad separada y enorme:
- ¿Dónde se guardan los FASTQs originales? ¿Están respaldados?
- ¿Qué ocurre con el directorio `work/` cuando el análisis termina? (Puede ocupar 10× más espacio que los resultados)
- `nextflow clean -f` — cuándo es seguro limpiar
- La diferencia entre datos de trabajo (work/), resultados finales (results/) y datos de archivo (long-term storage)

Sugeriría al menos 5 minutos en data management. Es algo que ningún predoc tiene organizado y que les ahorrará disgustos enormes.

#### El problema del "funciona en mi portátil pero no en el cluster"

Aunque Nextflow promete portabilidad, en la práctica hay problemas comunes que los estudiantes encontrarán:
- Diferentes versiones de Java en el cluster
- Singularity no disponible o con versión antigua
- Restricciones de acceso a internet desde los nodos de cómputo (no pueden descargar imágenes Docker)
- Límites de storage en el directorio home vs. scratch

No para desmotivar, sino para preparar. "Aquí están los 3 problemas más comunes y cómo resolverlos" es mucho más valioso que descubrirlos a las 11 de la noche antes de un deadline.

---

### Ajuste de timing: mi estimación es más realista que la tuya

Tus bloques suman 3-4 horas en papel, pero en la práctica:

- Las preguntas y discusiones siempre toman el doble de lo esperado.
- El setup técnico del hands-on para 20 personas es caótico sin preparación específica.
- La parte de "tus primeros pasos con Nextflow" no puede hacerse en 20 minutos si el objetivo es que realmente lo entiendan.

**Mi recomendación:**
- Si tienes 3 horas: corta el bloque 2 (comparativa) a 15 minutos y centra el tiempo en el bloque 4 (hands-on).
- Si tienes 4-5 horas: mantén la estructura pero añade los bloques de contenedores y ciclo de vida de un análisis real.
- Prepara un documento de "recursos para después" que cubra lo que no puedas dar en clase.

---

### Consejo final sobre la pedagogía del seminario

He visto muchos seminarios de bioinformática que fallan porque intentan transmitir conocimiento en lugar de **cambiar el comportamiento** de los asistentes. La pregunta correcta no es "¿qué he enseñado?", sino:

> **¿Qué va a hacer diferente un asistente a este seminario el lunes por la mañana?**

Si la respuesta es "va a instalar Nextflow y ejecutar `nextflow run nf-core/rnaseq -profile test,docker`", el seminario ha sido un éxito. Si la respuesta es "le parece interesante Nextflow y va a pensarlo", no.

Para ello, termina el seminario con un **call to action muy concreto**:
- Una lista de 3 pasos específicos para esta semana
- Un problema que tengan ahora mismo en su análisis y cómo Nextflow lo resolvería
- El compromiso de compartir en el Slack del curso en 2 semanas: "He ejecutado X con Nextflow"

---

## Resumen de recomendaciones del catedrático

| Recomendación | Prioridad | Impacto esperado |
|---|---|---|
| Añadir escenarios concretos ("el revisor 2 llegó") | Alta | Mayor conexión emocional con la audiencia |
| Inversión del orden: nf-core antes que Hello World | Alta | Reduce frustración, muestra el destino primero |
| Sección de errores comunes al empezar | Alta | Evita las frustraciones más comunes |
| Bloque de contenedores con más profundidad | Alta | Sin contenedores, la reproducibilidad es incompleta |
| 5 minutos sobre "cómo venderlo a tu PI" | Media | Elimina una barrera de adopción real |
| Abordar el miedo explícitamente | Media | Reduce la resistencia al cambio |
| Añadir bloque de gestión de datos | Media | Gap importante que nadie aborda |
| Ejemplo de ciclo de vida completo de análisis | Media | Conecta el Hello World con la práctica real |
| Control de versiones y trazabilidad | Media | Crucial para publicación científica |
| Ajustar el timing con márgenes realistas | Baja | Evita terminar con prisas el hands-on |

---

## Una última reflexión

Hay una asimetría que vale la pena nombrar explícitamente en clase: el valor de Nextflow es **asimétrico en el tiempo**.

- **Coste**: alto al principio (curva de aprendizaje, configuración inicial, cambio de hábitos)
- **Beneficio**: bajo al principio, **exponencialmente creciente** con el tiempo y la escala

El predoc que está procesando 10 muestras en su primer año no siente urgentemente la necesidad de Nextflow. El mismo investigador, 2 años después, procesando 200 muestras de scRNA-seq en colaboración con 3 grupos internacionales y 2 revisiones de paper pendientes, lo entiende perfectamente.

Tu seminario puede **acortar ese ciclo**: hacer que los estudiantes vean ese futuro hoy y empiecen a construir las habilidades antes de necesitarlas desesperadamente.

*Eso es exactamente lo que hace un buen profesor.*
