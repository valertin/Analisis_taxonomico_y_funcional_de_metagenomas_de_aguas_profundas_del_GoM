# Análisis taxonómico y funcional de metagenomas de aguas profundas del GoM

### Posdoctorante:
Dr. [Valentín Pérez-Hernández](https://www.researchgate.net/profile/Valentin-Perez-Hernandez-2)

[Laboratorio de Metagenómica](https://www.metagenomics-cicese.net/grupo)

<img width="3296" height="1546" alt="figure2-metagenomic" src="https://github.com/user-attachments/assets/a388c829-dbc8-4edc-923f-c01054823adc" />


## Departamento de Innovación Biomédica
## [Centro de Investigación Científica y de Educación Superior de Ensenada](https://www.cicese.edu.mx/) (CICESE)


## ``Introducción``

Se describe el proceso empleado para el análisis de datos metagenomicos provenientes de las aguas profundad del golfo de México. El procesamiento de las secuencias fueron procesadas siguiendo el flujo de trabajo indicado en la figura 1. 

<img width="4437" height="2538" alt="METAGENOMIC-Workflow" src="https://github.com/user-attachments/assets/7cf88fce-11fa-4b21-a6a3-1aa4e19dd558" />

Fig. 1. Flujo de trabajo del proceso de análisis bioinformático aplicado a las secuencias metagenomicas

# ``Pre-procesamiento``

El preprocesamiento en metagenómica es la etapa crítica de control de calidad y limpieza de las secuencias crudas (reads) obtenidas directamente del secuenciador, donde se eliminan adaptadores, primers, secuencias de baja calidad (con puntuaciones Phred bajas) y contaminantes externos, como ADN huésped (ej., humano o del organismo anfitrión). Este proceso utiliza herramientas como FastQC para diagnóstico, Sickle para recorte de reads, y BBMap o Bowtie2 para filtrar contaminantes mediante mapeo contra genomas de referencia.

Nuestras secuencias provienen de la plataforma Illumina Hi-Seq 2000, 2×100 paired end. Por lo que el procedimiento empleado indicará el empleo de los comandos para secuencias pareadas (paired end).

# Eliminacion de lecturas humanas

Las lecturas crudas (*R1 y *R2) se alinearon con [bbmap.sh](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) contra el hg19, de acuerdo con los parámetros publicados en Bushnell (2014).

## forward
```
$ bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/PATH/ref/ qtrim=rl trimq=10 untrim
-Xmx23g in=/PATH/R1.fastq outu=/PATH R1_cleaned.fastq outm=/R1_humanreads.fastq
```
## reverse
```
$ bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/PATH/ref/ qtrim=rl trimq=10 untrim
-Xmx23g in=/PATH/R2.fastq outu=/PATH R2_cleaned.fastq outm=/R2_humanreads.fastq
```
# Filtrado de calidad
El programa [Sickle](https://github.com/najoshi/sickle) se empleó para el control de calidad, truncando la calidad Phred a un valor de 20 (-q 20 es suficiente calidad).

```
$ sickle pe -f R1_cleaned.fastq -r R2_cleaned.fastq -t sanger -o R1_cleaned_filtered.fastq
-p R2_cleaned_filtered.fastq -s uneven_cleaned_filtered.fastq -n -q 20 -l 50
```
# ``Procesamiento``
# ``Anotación funcional``

## Ensamble: Megahit
El ensamble de los archivos se realizó empleando [MEGAHIT](https://github.com/voutcn/megahit); este permite emplear las secuencias que pasaron el control de calidad, pero que perdieron una lectura paired end (opción -r: uneven_cleaned_filtered.fastq). 

```
$ megahit --k-min 27 --k-max 101 --k-step 10 --min-contig-len 500 -1 R1_cleaned_filtered.fastq -2 R2_cleaned_filtered.fastq
-r uneven_cleaned_filtered.fastq -t $(nproc) -o assembly/
```

## Calculo de la cobertura de los contigs [BBMap](https://github.com/BioInfoTools/BBMap)
Para calcular la cobertura de los contigs con BBMap, se utiliza el principio de mapear las lecturas originales (reads) de vuelta a su propio ensamble. El comando bbmap.sh alinea los archivos FASTQ de lecturas (limpios o preprocesados) contra el archivo de contigs (en formato FASTA). La herramienta genera entonces un archivo BAM, que contiene la información de alineamiento para cada lectura. Finalmente, se usa el comando pileup.sh sobre el archivo BAM para generar un reporte de cobertura, que incluye la cobertura promedio por contig (la profundidad de secuenciamiento en cada posición) y la cobertura de bases (el porcentaje del contig cubierto por al menos una lectura).

```
$ bowtie2-build –threads 12 ../control/final.contigs.newheader.fa ctrl
$ samtools sort maiz.aligned.sam.bam > maiz.aligned.sam.bam.sorted.bam
$ samtools index control.aligned.sam.bam.sorted.bam
```

## Calcula información de cobertura por contigs [pileup.sh](https://nf-co.re/modules/bbmap_pileup/) 
Después de mapear las lecturas contra tus contigs con bbmap.sh (que genera un archivo BAM), ejecutas el comando pileup.sh de la suite BBMap usando ese archivo BAM como entrada. Este comando procesa el archivo de alineamiento y calcula automáticamente métricas clave por contig, como la cobertura promedio (Avg_fold), el porcentaje de bases cubiertas (Covered_percent), y la longitud, generando un reporte tabular listo para su análisis posterior.

```
$ pileup.sh in=mapping/aligned.sam.gz out=mapping/coverage.txt 2>&1  | tee mapping/log_coverage.log
```

## Calculo de metricas de los ensambles (contigs) con [MetaQuast](https://quast.sourceforge.net/metaquast)
MetaQuast evalúa la calidad de los ensambles generados por MEGAHIT (u otros ensambladores), calculando métricas clave como el tamaño total del ensamble, el número de contigs, la longitud de N50 y la detección de errores como misassemblies. Proporciona reportes HTML interactivos y tablas resumen que permiten visualizar y comparar datos estadisticos del ensamble.

```
metaquast.py -o metaquast/final.contigs.fa
```

## Anotación funcional con [PROKKA](https://github.com/tseemann/prokka)
La anotación funcional con PROKKA es un proceso automatizado que predice y describe genes en secuencias de contigs, identificando elementos genómicos como genes codificantes de proteínas (CDS), ARN ribosomal (rRNA), ARN de transferencia (tRNA) y otros elementos funcionales. Utiliza la herramienta Prodigal (para predicción de genes) y bases de datos como UniProt, Pfam y NCBI para asignar nombres de genes, productos génicos y números EC (Enzyme Commission), generando archivos de salida estandarizados (GBK, GFF, FAA)

```
$ Prokka final.contigs.fa  --outdir treatment/ --prefix control_ --metagenome --norrna --notrna
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
kraken2 --db $KRAKEN_DB --paired --threads 24 --report $OUTPUT_DIR/sample_report.txt --output $OUTPUT_DIR/sample_output.kraken
  $READS_DIR/sample_R1.fastq.gz   $READS_DIR/sample_R2.fastq.gz
```

# Re-estimación de la abundancia con [Bracken](https://github.com/jenniferlu717/Bracken)
Bracken (Bayesian Reestimation of Abundance with KrakEN) toma el reporte de Kraken2  y aplica un modelo bayesiano para recalcular y refinar las abundancias relativas de cada taxón.

```
bracken -d $KRAKEN_DB -i $OUTPUT_DIR/sample_report.txt -o $OUTPUT_DIR/sample_bracken.species -l S -t 50
```

# Análisis estadisticos y visualización con [R](https://www.r-project.org/)

```
Uso de R y R studio
```

## Conclusiones

El flujo de trabajo presentado abarca desde el control de calidad de las lecturas hasta la clasificación taxonómica (con herramientas como Kraken2 y Bracken), ensamblaje (con MEGAHIT), anotación funcional (PROKKA, BlastKOALA), se sugiere el uso de herramientas como R para análisis diferencial (DESeq2) y visualización (redes, heatmaps, PCA/PCoA). Los códigos proporcionados son ejemplos generales que requieren ajustes específicos según las características de cada dataset (calidad de lecturas, complejidad de la comunidad, objetivos del estudio). Este pipeline sirve como marco de referencia, pero su aplicación exitosa depende de la optimización para cada caso particular y la integración de evidencias multidisciplinarias.

### Referencias

- Andrews, S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data. 
- Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, August 2014, Pages 2114–2120, https://doi.org/10.1093/bioinformatics/btu170
- Dinghua Li, Chi-Man Liu, Ruibang Luo, Kunihiko Sadakane, Tak-Wah Lam, MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph, Bioinformatics, Volume 31, Issue 10, May 2015, Pages 1674–1676, https://doi.org/10.1093/bioinformatics/btv033
- Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104 https://doi.org/10.7717/peerj-cs.104
- R Core Team (2025). _R: A Language and Environment for Statistical computing
- Torsten Seemann, Prokka: rapid prokaryotic genome annotation, Bioinformatics, Volume 30, Issue 14, July 2014, Pages 2068–2069, https://doi.org/10.1093/bioinformatics/btu153
- Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
