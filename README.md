# Análisis taxonómico y funcional de metagenomas de aguas profundas del GoM

### Posdoctorante:
Dr. [Valentín Pérez-Hernández](https://www.researchgate.net/profile/Valentin-Perez-Hernandez-2)

[Laboratorio de Metagenómica](https://www.metagenomics-cicese.net/grupo)

<img width="3296" height="1546" alt="figure2-metagenomic" src="https://github.com/user-attachments/assets/a388c829-dbc8-4edc-923f-c01054823adc" />


## Departamento de Innovación Biomédica
## [Centro de Investigación Científica y de Educación Superior de Ensenada](https://www.cicese.edu.mx/) (CICESE)


## ``Introducción``

El estudio de las comunidades microbianas en ambientes extremos, como las aguas profundas, representa un desafío técnico y biológico debido a las condiciones de alta presión, baja temperatura y escasez de nutrientes que caracterizan estos ecosistemas. La metagenómica ha emergido como una herramienta fundamental para explorar la diversidad y funcionalidad de microorganismos no cultivables, permitiendo descifrar su potencial biotecnológico y ecológico. En este trabajo, se describe el proceso empleado para el análisis de datos metagenómicos provenientes de muestras de aguas profundas, utilizando un flujo de trabajo computacional que abarca desde el control de calidad de las secuencias hasta la anotación taxonómica y funcional. Mediante el uso de herramientas como Kraken2 para la clasificación taxonómica, MEGAHIT para el ensamblaje de novo, y PROKKA junto con KEGG para la predicción de funciones metabólicas, se busca caracterizar la composición microbiana y su papel en los ciclos biogeoquímicos en este entorno único. 

En este documento se describe de forma general el proceso empleado para el análisis de datos metagenomicos provenientes de las aguas profundas del golfo de México. El procesamiento de las secuencias fueron realizadas siguiendo el flujo de trabajo indicado en la figura 1. 

<img width="4437" height="2538" alt="METAGENOMIC-Workflow" src="https://github.com/user-attachments/assets/7cf88fce-11fa-4b21-a6a3-1aa4e19dd558" />

Fig. 1. Flujo de trabajo del proceso de análisis bioinformático aplicado a las secuencias metagenomicas

# ``Pre-procesamiento``

El preprocesamiento en metagenómica es la etapa crítica de control de calidad y limpieza de las secuencias crudas (reads) obtenidas directamente del secuenciador, donde se eliminan adaptadores, primers, secuencias de baja calidad (con puntuaciones Phred bajas) y contaminantes externos, como ADN huésped (ej., humano o del organismo anfitrión). Este proceso utiliza herramientas como FastQC para diagnóstico, Sickle para recorte de reads, y BBMap o Bowtie2 para filtrar contaminantes mediante mapeo contra genomas de referencia.

Nuestras secuencias provienen de la plataforma Illumina Hi-Seq 2000, 2×100 paired end. Por lo que el procedimiento empleado indicará el empleo de los comandos para secuencias pareadas (paired end).

# Eliminacion de lecturas humanas

Las lecturas crudas (*R1 y *R2) se alinearon con [bbmap.sh](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) contra el hg19, de acuerdo con los parámetros publicados en Bushnell (2014).

## Comando para forward y reverse
```
$ bbmap.sh \
  minid=0.95
  maxindel=3
  bwr=0.16
  bw=12
  quickmatch
  fast
  minhits=2
  path=/PATH/ref/
  qtrim=rl
  trimq=10
  untrim
  -Xmx23g
  in=/PATH/R1.fastq
  outu=/PATH R1_cleaned.fastq
  outm=/R1_humanreads.fastq
```

```
Detalles del comando:

* minid:	Establece la identidad mínima de alineamiento en 95%. Solo se consideran lecturas que se alinean con al menos 95% de similitud.
* maxindel:	Permite un máximo de 3 inserciones/deleciones por alineamiento. Mejora la precisión en regiones conservadas.
* bwr: Define el ratio de retroceso en el algoritmo de alineamiento. Valores bajos hacen el mapeo más estricto.
* bw:	Establece el ancho de banda para el alineamiento. Controla cuántas bases se consideran en cada paso.
* quickmatch:	Activa un modo rápido que acelera el mapeo usando coincidencias simplificadas.
* fast:	Usa una heurística rápida para acelerar el proceso, sacrificando algo de sensibilidad.
* minhits:	Requiere al menos 2 k-mers coincidentes para considerar una lectura como mapeada.
* path:	Ruta a la carpeta que contiene el índice de referencia (por ejemplo, genoma humano).
* qtrim:	Recorta las lecturas en ambos extremos (right y left) si la calidad cae por debajo del umbral.
* trimq:	Umbral de calidad para el recorte. Bases con calidad <10 serán eliminadas en los extremos.
* untrim:	Indica que las lecturas recortadas deben ser reconstruidas a su forma original tras el mapeo.
* -Xmx:	Asigna 23 GB de memoria RAM para el proceso (ajusta según tu sistema).
* in: Archivo de entrada con las lecturas crudas.
* outu:	Archivo de salida con las lecturas no mapeadas.
* outm:	Archivo de salida con las lecturas mapeadas.
```
# Filtrado de calidad
El programa [Sickle](https://github.com/najoshi/sickle) se empleó para el control de calidad, truncando la calidad Phred a un valor de 20 (-q 20 es suficiente calidad).

```
$ sickle pe \
  -f R1_cleaned.fastq \
  -r R2_cleaned.fastq \
  -t sanger \
  -o R1_cleaned_filtered.fastq \
  -p R2_cleaned_filtered.fastq \
  -s uneven_cleaned_filtered.fastq \
  -n \
  -q 20 \
  -l 50
```

```
Detalles de los parámetro:
* pe:	Modo paired-end: procesa lecturas emparejadas.
* -f: Archivo de lectura forward (R1) ya limpiado.
* -r: Archivo de lectura reverse (R2) ya limpiado.
* -t: Formato de calidad de las lecturas (Phred+33, típico en Illumina).
* -o: Salida de lecturas forward que pasaron el filtro.
* -p: Salida de lecturas reverse que pasaron el filtro.
* -s: Lecturas que quedaron sin pareja tras el filtrado (por ejemplo, si R1 pasa pero R2 no).
* -n: No convierte las bases de baja calidad en "N"; simplemente las recorta.
* -q: Umbral de calidad, recorta bases con calidad <20 (Phred score).
* -l: Longitud mínima, descarta lecturas que queden con menos de 50 bases tras el recorte.
```

# ``Procesamiento``
# ``Anotación funcional``

## Ensamble: Megahit
El ensamble de los archivos se realizó empleando [MEGAHIT](https://github.com/voutcn/megahit); este permite emplear las secuencias que pasaron el control de calidad, pero que perdieron una lectura paired end. 

```
$ megahit \
  --k-min 27 \
  --k-max 101 \
  --k-step 10 \
  --min-contig-len 500 \
  -1 R1_cleaned_filtered.fastq \
  -2 R2_cleaned_filtered.fastq \
  -r uneven_cleaned_filtered.fastq \
  -t $(nproc) \
  -o assembly/
```

```
Detalles de los parámetros
--k-min:	Tamaño mínimo de k-mer usado en el ensamblaje.
--k-max:	Tamaño máximo de k-mer.
--k-step:	Incremento entre tamaños de k-mer (usa 27, 37, 47... hasta 101).
--min-contig-len: Filtra contigs menores a 500 bp en la salida final.
-1: Lecturas forward (paired-end) ya filtradas.
-2: Lecturas reverse (paired-end) ya filtradas.
-r: Lecturas single-end (sin pareja tras filtrado con Sickle).
-t: Usa todos los núcleos disponibles del procesador para acelerar el ensamblaje.
-o: Carpeta de salida donde se guardan los contigs ensamblados.
```

## Calculo de la cobertura de los contigs [BBMap](https://github.com/BioInfoTools/BBMap)
Para calcular la cobertura de los contigs con BBMap, se utiliza el principio de mapear las lecturas originales (reads) de vuelta a su propio ensamble. El comando bbmap.sh alinea los archivos FASTQ de lecturas (limpios o preprocesados) contra el archivo de contigs (en formato FASTA). La herramienta genera entonces un archivo BAM, que contiene la información de alineamiento para cada lectura. Finalmente, se usa el comando pileup.sh sobre el archivo BAM para generar un reporte de cobertura, que incluye la cobertura promedio por contig (la profundidad de secuenciamiento en cada posición) y la cobertura de bases (el porcentaje del contig cubierto por al menos una lectura).

```
$ bowtie2-build –threads 12 ../control/final.contigs.newheader.fa ctrl
```
```
Detalles de los parametros:
* bowtie2-build:	Ejecuta el constructor de índice de Bowtie2.
* --threads: Núcleos de CPU usados en el proceso.
* final.contigs.fa:	Archivo FASTA que contiene los contigs ensamblados. Este es el archivo que se convertirá en índice.
* ctrl:	Prefijo del índice de salida. Bowtie2 generará varios archivos con nombres como ctrl.1.bt2, ctrl.2.bt2, etc., que conforman el índice completo.
```
```
$ samtools sort sample.aligned.sam.bam > sample.aligned.sam.bam.sorted.bam
```
```
Detalles de los parametros:
* samtools sort:	Ordena el archivo BAM por posición en el genoma de referencia.
* sample.aligned.sam.bam:	Archivo BAM de entrada, generado previamente por Bowtie2 u otra herramienta de mapeo.
* >:	Redirecciona la salida estándar al archivo especificado.
* sample.aligned.sam.bam.sorted.bam:	Archivo BAM de salida, ya ordenado por coordenadas genómicas.
```

```
$ samtools index control.aligned.sam.bam.sorted.bam -o output
```

```
Detalles de los parametros:
* index	Subcomando que genera un archivo de índice .bai para un archivo BAM.
* control.aligned.sam.bam.sorted.bam	Archivo BAM de entrada, que debe estar ordenado por coordenadas (output típico de samtools sort).
* -o:	Permite especificar el nombre del archivo de índice manualment
```
## Calcula información de cobertura por contigs [pileup.sh](https://nf-co.re/modules/bbmap_pileup/) 
Después de mapear las lecturas contra tus contigs con bbmap.sh (que genera un archivo BAM), ejecutas el comando pileup.sh de la suite BBMap usando ese archivo BAM como entrada. Este comando procesa el archivo de alineamiento y calcula automáticamente métricas clave por contig, como la cobertura promedio (Avg_fold), el porcentaje de bases cubiertas (Covered_percent), y la longitud, generando un reporte tabular listo para su análisis posterior.

```
$ pileup.sh in=mapping/aligned.sam.gz out=mapping/coverage.txt 
```

```
Detalles de los parametros:
* in=	Archivo de entrada con las lecturas alineadas (en formato SAM comprimido). Este archivo debe haber sido generado previamente por bbmap.sh o una herramienta de mapeo.
* out=	Archivo de salida que contendrá las estadísticas de cobertura por contig. Se genera en formato texto tabulado.
```
## Calculo de metricas de los ensambles (contigs) con [MetaQuast](https://quast.sourceforge.net/metaquast)
MetaQuast evalúa la calidad de los ensambles generados por MEGAHIT (u otros ensambladores), calculando métricas clave como el tamaño total del ensamble, el número de contigs, la longitud de N50 y la detección de errores como misassemblies. Proporciona reportes HTML interactivos y tablas resumen que permiten visualizar y comparar datos estadisticos del ensamble.

```
metaquast.py final.contigs.fa -o metaquast/
```

```
Detalles de los parametros:
* final.contigs.fa: secuencias ensambladas
* -o: archivos de salida con los estadisticos de los contigs evaluados
```
## Anotación funcional con [PROKKA](https://github.com/tseemann/prokka)
La anotación funcional con PROKKA es un proceso automatizado que predice y describe genes en secuencias de contigs, identificando elementos genómicos como genes codificantes de proteínas (CDS), ARN ribosomal (rRNA), ARN de transferencia (tRNA) y otros elementos funcionales. Utiliza la herramienta Prodigal (para predicción de genes) y bases de datos como UniProt, Pfam y NCBI para asignar nombres de genes, productos génicos y números EC (Enzyme Commission), generando archivos de salida estandarizados (GBK, GFF, FAA)

```
$ prokka final.contigs.fa \
  --outdir treatment/ \
  --prefix control_ \
  --metagenome \
  --norrna \
  --notrna
```

```
Detalles de los parametros:
* final.contigs.fa	Archivo FASTA con los contigs ensamblados que se van a anotar.
* --outdir Directorio donde se guardarán los archivos de salida.
* --prefix Prefijo para nombrar los archivos generados (ej. control_.gff, control_.faa, etc.).
* --metagenome	Indica que el ensamblaje proviene de un metagenoma, lo que ajusta el algoritmo para tratar contigs fragmentados y sin orden definido.
* --norrna	Omite la búsqueda de genes de ARN ribosomal (rRNA). Útil si ya se han filtrado o si no se desea incluirlos.
* --notrna	Omite la búsqueda de genes de ARN de transferencia (tRNA).
```

## Análisis de vías metabólicas [KofamKOALA-BlastKoala](https://www.genome.jp/tools/kofamkoala/) 
KofamKOALA  toma la anotación funcional de PROKKA (específicamente el archivo FAA de secuencias proteicas predichas) y realiza búsquedas de similitud de secuencias contra la base de datos KEGG, asignando números K y reconstruyendo pathways metabólicos. Al analizar los genes anotados por PROKKA, BlastKOALA contextualiza funcionalmente los contigs al mapear enzimas y reacciones bioquímicas en redes de KEGG, permitiendo interpretar roles metabólicos, ciclos biogeoquímicos o potencial biotecnológico 

```
Web base
```
<img width="795" height="743" alt="image" src="https://github.com/user-attachments/assets/74e05e6c-90de-4742-b754-3580fd667ac7" />

# Anotación taxónomica con [Kraken2](https://github.com/DerrickWood/kraken2)
Kraken2 asigna taxonomía a los reads de metagenoma mediante un método de coincidencia exacta de k-mers. Para cada lectura, divide su secuencia en pequeños fragmentos de longitud fija (k-mers, por ejemplo k=35) y consulta cada k-mer en una base de datos preconstruida que contiene genomas de referencia asociados a taxones específicos.

```
$ kraken2
* --db $KRAKEN_DB
* --paired
* --threads 24
* --report $OUTPUT_DIR/sample_report.txt
* --output $OUTPUT_DIR/sample_output.kraken
* $READS_DIR/sample_R1.fastq.gz
* $READS_DIR/sample_R2.fastq.gz
```

```
* --db Ruta a la base de datos de Kraken2 (puede ser MiniKraken, RefSeq, etc.)
* --paired	Indica que las lecturas son pareadas (R1 y R2)
* --threads	nucleos a emplear durante el procesamiento de los datos
* --report Genera un resumen tabular con abundancias por taxón
* --outputArchivo con la clasificación individual de cada lectura
* $READS_DIR/sample_R1.fastq.gz Archivos FASTQ de lecturas R1
* $READS_DIR/sample_R2.fastq.gz	Archivos FASTQ de lecturas R2
```
# Re-estimación de la abundancia con [Bracken](https://github.com/jenniferlu717/Bracken)
Bracken (Bayesian Reestimation of Abundance with KrakEN) toma el reporte de Kraken2  y aplica un modelo bayesiano para recalcular y refinar las abundancias relativas de cada taxón.

```
$ bracken \
  -d $KRAKEN_DB \
  -i $OUTPUT_DIR/sample_report.txt \
  -o $OUTPUT_DIR/sample_bracken.species \
  -l S \
  -t 50
```

```
* -d: Ruta a la base de datos de Kraken2 usada para la clasificación. Bracken necesita esta base para conocer la distribución de k-mers por taxón.
* -i: rchivo de entrada generado por Kraken2 con conteos por taxón.
* -o:	Archivo de salida con las abundancias corregidas por Bracken.
* -l: Nivel taxonómico para el análisis. S indica especie (puedes usar G para género, F para familia, etc.).
* -t: Número máximo de lecturas consideradas por taxón para el reestimado. Ayuda a controlar el tiempo de ejecución y la precisión.
```

# Análisis estadisticos y visualización con [R](https://www.r-project.org/)

```
Uso de R y R studio
```
![protozoa+pie](https://github.com/user-attachments/assets/71aaae5b-5cd6-4756-8983-9bf5fcbd66cb)  <img width="740" height="605" alt="KEGG2+PCA" src="https://github.com/user-attachments/assets/d3651cc6-4d39-4263-8372-6634141255c7" />

Fig. 2 Ejemplos de gráficos que se pueden realizar en R con los resultados obtenidos del análisis de secuencias.

## Conclusiones

El flujo de trabajo presentado abarca desde el control de calidad de las lecturas hasta la clasificación taxonómica (con herramientas como Kraken2 y Bracken), ensamblaje (con MEGAHIT), anotación funcional (PROKKA, BlastKOALA), se sugiere el uso de herramientas como R para análisis diferencial (DESeq2) y visualización (redes, heatmaps, PCA/PCoA). Los códigos proporcionados son ejemplos generales que requieren ajustes específicos según las características de cada dataset (calidad de lecturas, complejidad de la comunidad, objetivos del estudio). Este pipeline sirve como marco de referencia, pero su aplicación exitosa depende de la optimización para cada caso particular.

### Referencias

- Andrews, S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data. 
- Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, August 2014, Pages 2114–2120, https://doi.org/10.1093/bioinformatics/btu170
- Dinghua Li, Chi-Man Liu, Ruibang Luo, Kunihiko Sadakane, Tak-Wah Lam, MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph, Bioinformatics, Volume 31, Issue 10, May 2015, Pages 1674–1676, https://doi.org/10.1093/bioinformatics/btv033
- Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104 https://doi.org/10.7717/peerj-cs.104
- R Core Team (2025). _R: A Language and Environment for Statistical computing
- Torsten Seemann, Prokka: rapid prokaryotic genome annotation, Bioinformatics, Volume 30, Issue 14, July 2014, Pages 2068–2069, https://doi.org/10.1093/bioinformatics/btu153
- Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
