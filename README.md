# Análisis taxonómico y funcional de metagenomas de aguas profundas del GoM

### Posdoctorante:
Dr. [Valentín Pérez-Hernández](https://www.researchgate.net/profile/Valentin-Perez-Hernandez-2)

[Laboratorio de Metagenómica](https://www.metagenomics-cicese.net/grupo)

<img width="3296" height="1546" alt="figure2-metagenomic" src="https://github.com/user-attachments/assets/a388c829-dbc8-4edc-923f-c01054823adc" />


## Departamento de Innovación Biomédica
## [Centro de Investigación Científica y de Educación Superior de Ensenada](https://www.cicese.edu.mx/) (CICESE)


## ``Introducción``

Se describe el proceso empleado para el análisis de datos metagenomicos provenientes de las aguas profundad del golfo de México. El procesamiento de las secuencias fueron procesadas siguiendo el flujo de trabajo indicado en la figura 1. 

<img width="4437" height="2538" alt="METAGENOMIC-Workflow" src="https://github.com/user-attachments/assets/6438f462-402c-4656-8906-af46224e5a2e" />
Fig. 1. Flujo de trabajo del proceso de análisis bioinformático aplicado a las secuencias metagenomicas

# ``Pre-procesamiento``

El preprocesamiento en metagenómica es la etapa crítica de control de calidad y limpieza de las secuencias crudas (reads) obtenidas directamente del secuenciador, donde se eliminan adaptadores, primers, secuencias de baja calidad (con puntuaciones Phred bajas) y contaminantes externos, como ADN huésped (ej., humano o del organismo anfitrión). Este proceso utiliza herramientas como FastQC para diagnóstico, Sickle para recorte de reads, y BBMap o Bowtie2 para filtrar contaminantes mediante mapeo contra genomas de referencia.

Nuestras secuencias provienen de la plataforma Illumina Hi-Seq 2000, 2×100 paired end. Por lo que el procedimiento empleado indicará el empleo de los comandos para secuencias pareadas (paired end).

# Eliminacion de lecturas humanas

Las lecturas crudas (*R1 y *R2) se alinearon con bbmap.sh contra el hg19, de acuerdo con los parámetros publicados en Bushnell (2014).

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
$ sickle pe -f R1_cleaned.fastq -r R2_cleaned.fastq -t sanger -o R1_cleaned_filtered.fastq -p R2_cleaned_filtered.fastq -s uneven_cleaned_filtered.fastq -n -q 20 -l 50
```
