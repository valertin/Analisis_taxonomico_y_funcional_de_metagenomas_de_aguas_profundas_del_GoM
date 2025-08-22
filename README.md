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

# Anotación taxónomica

