# Biol217 Practice Session Day-1

  what we have earned so far?

    1. Basic Linux
    2. Bioinformatics basic understandig
    3. Linux commands

    -   copy from one folder to another_
  

    Block of code:

  ```sh
  cp source destination
  ```

  inline code:
  this is the command `cp`

damit ich ins caucluster komme
  ssh -X sunam238@caucluster.rz.uni-kiel.de

### Task: How to add links and images to the end file

-------

# Day 2

##  code for Sequencing files

zuerst muss gesagt werden, was der Umstand für das Programmieren ist. Wir haben in dem Fall jetzt sequence files, 6 Stück insgesamt, 3 davon sind normal und 3 davon reversed reads. Am Ende ergeben 1 normal read und ein reversed read paired end readfiles. Die Files habe den Datentypen fastq.gz.
Mit fastqc (qc für QualityControl) verwirklicht man eine Visualisierung der fastq Dateien. HTMLs werden ausgegeben, die man dann in einen Internetbrowser kopieren kann und die visualisierte Datei sieht. 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/0_raw_reads


for i in *.gz; do fastqc $i -o ../1_fastqc/; done
```
--------------------

```
dann muss man gucken, ob die File vom caucluster prozessiert wurde: 

squeue -u sunam238
ls 
queue 
```
Die html s auf den eigenen PC kopieren, dafür einfach neues Terminal-Fenster öffnen, denn das aktuelle Terminal greift auf das cau-cluster zu.
```

scp sunam238@caucluster.rz.uni-kiel.de:/work_beegfs/suname238/Metagenomics/1_fastqc/*.html.
```
Mit fastp cleaned man die reads und gibt neuen input

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/0_raw_reads
```

## fastp 

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8


fastp -i BGR_130305_mapped_R1.fastq.gz  -I BGR_130305_mapped_R2.fastq.gz  -R BGR_130305_reportfile -o ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130305_mapped_R2_clean.fastq.gz -t 6 -q 20

fastp -i BGR_130527_mapped_R1.fastq.gz  -I BGR_130527_mapped_R2.fastq.gz  -R BGR_130527_reportfile -o ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -t 6 -q 20

fastp -i BGR_130708_mapped_R1.fastq.gz  -I BGR_130708_mapped_R2.fastq.gz  -R BGR_130708_reportfile -o ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -t 6 -q 20

```
---------------------

## megahit

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=megahit.out
#SBATCH --error=megahit.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/2_fastp
                                       
megahit -1 GR_130305_mapped_R1.fastq.gz -1 BGR_130527_mapped_R1.fastq.gz -1 BGR_130708_mapped_R1.fastq.gz -2 GR_130305_mapped_R2.fastq.gz -2 BGR_130527_mapped_R2.fastq.gz -2 BGR_130708_mapped_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o ../3_coassembly -t 12 
```
Day 3

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics

metaquast -t 6 -o ./3_metaquast -m 1000 ./3_coassembly/final.contigs.fa

----
danach ein zweites Fenster im Terminal öffnen und mit dem stu-Account den Befehl eingeben zum Kopieren (scp=):
----
scp sunam238@caucluster.rz.uni-kiel.de:/work_beegfs/sunam238/Metagenomics/3_metaquast/report.html .



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8


anvi-script-reformat-fasta ./3_coassembly/final.contigs.fa -o /work_beegfs/sunam238/Metagenomics/3_coassembly/contigs.anvio.fa --min-len 1000 --simplify-names --report-file name_conversion.txt

sbatch

squeue -u sunam238



----------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=bowtie
#SBATCH --output=bowtie.out
#SBATCH --error=bowtie.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8


cd work_beegfs/sunam238/Metagenomics/3_coassembly/

module load bowtie2
bowtie2-build contigs.anvio.fa contigs.anvio.fa.index

  
  
 


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -2 /PATH/TO/BGR_130305_mapped_R2_clean.fastq.gz -S map130305

bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -S map130527.sam

bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -S map130708.sam


-----
module load samtools
samtools view -bS ? > bam_file.bam

-------- 
Contigs data preparation
------

cd /work_beegfs/sunam238/Metagenomics/3_coassembly/


anvi-gen-contigs-database -f contigs.anvio.fa -o contigs.db -n 'biol217'


anvi-run-hmms -c contigs.db


srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --nodelist=node002 /bin/bash

---
------
----
for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done


anvi-profile -i ? -c ? --output-dir ?
----

Day 4

cd /workbeegfs/Metagenomics/3_coassembly

module load samtools

for i in *.sam; do samtools view -bS $i > "$i".bam; done

dieser Schritt war erforderlich, um aus den sam-Dateitypen (sequence mapping file) einen bam-Dateitypen (binary alignment file) zu machen.

-----

Fortsetzung von dem letzten Command von Day3.
Eine File wurde heruntergeladen -> contigs.db darin sind die contigs enthalten

```
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --nodelist=node002 /bin/bash 
```

Danach:

Muss ich auf anvio interactive zugreifen
anvio interactive brauch ich jedes einzelne Mal
ich muss immer die gleichen Schritte gehen
Ich muss die command line, die ich im interactive mode runnen möchte, replacen 
Das folgende Skript dient zur Visualisierung unserer Dateien.

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8
```
```
anvi-display-contigs-stats contigs.db
```
(Dieser command ist abhängig davon, was man im interactive mode runnen möchte und muss demnach geändert bleiben. Der Rest des Skriptes ist fix und kann für andere Visualisierungen verwendet werden.)

Neues Terminal öffnen!

ssh -L 8060:localhost:8080 sunam238@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node238

Internetbrowser öffnen und folgendes kopieren: 
http://127.0.0.1:8060 oder http://127.0.0.1:8080

anvi-display-contigs-stats contigs.db
(without srun, das klappt auch)

### Binning mit ANVI´O

ANVI´O: ANalysis and Visualization platform for microbial ´Omics

Online-Tutorial, wenn ich verwirrt bin: https://merenlab.org/2016/06/22/anvio-tutorial-v2/

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=bowtie
#SBATCH --output=bowtie.out
#SBATCH --error=bowtie.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8


cd /work_beegfs/sunam238/Metagenomics/3_coassembly/ -->




for i in *.bam; do anvi-init-bam $i -o ../5_anvio_profiles/"$i".sorted.bam; done
```



cd /work_beegfs/sunam238/Metagenomics/5_anvio_profiles

anvi-profile -i anvio.bam -c ../3_coassembly/contigs.db --output-dir OUTPUT_DIR

stattdessen das für jede Bam-File einzeln einzugeben, ist hier die For-Schleife:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=profiling
#SBATCH --output=profiling.out
#SBATCH --error=profiling.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/5_anvio_profiles

mkdir /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/profiling

for i in `ls *.sorted.bam | cut -d "." -f 1`; do anvi-profile -i "$i".bam.sorted.bam -c ../3_coassembly/contigs.db -o /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/profiling/”$i”; done
```

nächster Schritt: Verschmelzen der Anvi-Profile von verschiedenen Samples zu einem einzigen Profil:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=anvimerge
#SBATCH --output=anvimerge.out
#SBATCH --error=anvimerge.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/

anvi-merge ./5_anvio_profiles/profiling/map130305/PROFILE.db ./5_anvio_profiles/profiling/map130305/PROFILE.db ./5_anvio_profiles/profiling/map130708/PROFILE.db -o ./5_anvio_profiles/profiling/merged_profiles -c ./3_coassembly/contigs.db --enforce-hierarchical-clustering 
```
/5_anvio_profiles/profiling/map130305/
/5_anvio_profiles/profiling/map130527/
/5_anvio_profiles/profiling/map130708/


# danach kann mit dem Binning begonnen werden:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=anvicluster
#SBATCH --output=anvicluster.out
#SBATCH --error=anvicluster.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/

anvi-cluster-contigs -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -C METABAT --driver metabat2 --just-do-it --log-file log-metabat2

anvi-summarize -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -o SUMMARY_METABAT -C METABAT
```

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=maxbin
#SBATCH --output=maxbin.out
#SBATCH --error=maxbin.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/

anvi-cluster-contigs -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2

anvi-summarize -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -o MAXBIN2 -C MAXBIN2
```

### MAGs Quality Estimation

```
anvi-estimate-genome-completeness -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -C METABAT2
```

```
anvi-estimate-genome-completeness -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db --list-collections
```



```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-interactive -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -C YOUR_COLLECTION
```


# Day 4

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=maxbin
#SBATCH --output=maxbin.out
#SBATCH --error=maxbin.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics

anvi-summarize -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db --list-collections

anvi-summarize -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/profile.db -C METABAT2 -o SUMMARY_METABAT2 --just-do-it
```
anvi-estimate-genome-completeness -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -C METABAT > file

(damit lässt man sich eine Tabelle ausgeben, jetzt müssen die Archeae rausgepickt werden)

```
METABAT__23 | ARCHAEA
METABAT__44 | ARCHAEA
METABAT__1  | ARCHAEA 
```


```
cd /CC/5_anvio_profiles/MAXBIN2/bin_by_bin

mkdir ../../ARCHAEA_BIN_REFINEMENT

cp /PATH/TO/BIN_FOLDER_INFO_FROM_ERR_FILE/*.fa /PATH/TO/ARCHAEA_BIN_REFINEMENT/
```

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate gunc
```

```
cd /PATH/TO/ARCHAEA_BIN_REFINEMENT

mkdir GUNC

for i in *.fa; do gunc run -i ? -r /home/sunam238/Databases/gunc_db_progenomes2.1.dmnd --out_dir ? --threads 10 --detailed_output; done

```

anvi-summarize -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db --list-collections

anvi-summarize -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -C BINSANITY -o ./CC/5_anvio_profiles/BINSANITY --just-do-it > binsanity_file.sh



```
cd ./CC/5_anvio_profiles/SUMMARY_METABAT2/bin_by_bin

mkdir ../../ARCHAEA_BIN_REFINEMENT

cp METABAT_23/*.fa /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT/
cp METABAT_44/*.fa /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT/
cp METABAT_1/*.fa /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT/
```


# GUNC

```

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=gunc
#SBATCH --output=gunc.out
#SBATCH --error=gunc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate gunc

cd /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT


mkdir GUNC

gunc run -i METABAT__44-contigs.fa -r /work_beegfs/sunam238/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC/METABAT__44 --threads 10 --detailed_output

gunc run -i METABAT__23-contigs.fa -r /work_beegfs/sunam238/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC/METABAT__23 --threads 10 --detailed_output

gunc run -i METABAT__1-contigs.fa -r /work_beegfs/sunam238/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC/METABAT__1 --threads 10 --detailed_output

```


```
anvi-refine -c ./CC/4_mapping/contigs.db -C METABAT -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db --bin-id METABAT__44
```



gunc plot -d GUNC/METABAT__1/diamond_output/METABAT__1-contigs.diamond.progenomes_2.1.out -g GUNC/METABAT__1/gene_calls/gene_counts.json 





```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=day5
#SBATCH --output=day5.out
#SBATCH --error=day5.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT

anvi-run-scg-taxonomy -c ./CC/4_mapping/contigs.db -T 20 -P 2



```


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=hope5
#SBATCH --output=hope.out
#SBATCH --error=hope.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/

#CC/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT

anvi-estimate-scg-taxonomy -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy

anvi-estimate-scg-taxonomy -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt

anvi-summarize -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -o ./CC/5_anvio_profiles/merged_profiles/SUMMARY_METABAT2 -C METABAT



```



# Day6

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=01_fastqc
#SBATCH --output=01_fastqc.out
#SBATCH --error=01_fastqc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load micromamba/1.4.2
micromamba activate 01_short_reads_qc


# create new folder for output of qc 
mkdir -p $WORK/genomics/1_short_reads_qc/1_fastqc_raw
for i in $WORK/genomics/0_raw_reads_/short_reads/*.gz; do fastqc $i -o .$WORK/genomics/1_short_reads_qc/1_fastqc_raw -t 32; done

jobinfo
```


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=02_long_reads_qc
#SBATCH --output=02_long_reads_qc.out
#SBATCH --error=02_long_reads_qc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"

echo "---------long reads cleaning started---------"
micromamba activate 02_long_reads_qc

## 2.1 Nanoplot raw
cd $WORK/genomics/0_raw_reads/long_reads/
mkdir -p $WORK/genomics/2_long_reads_qc/1_nanoplot_raw
NanoPlot --fastq $WORK/genomics/0_raw_reads/long_reads/*.gz \
 -o $WORK/genomics/2_long_reads_qc/1_nanoplot_raw -t 32 \
 --maxlength 40000 --minlength 1000 --plots kde --format png \
 --N50 --dpi 300 --store --raw --tsv_stats --info_in_report

## 2.2 Filtlong
mkdir -p $WORK/genomics/2_long_reads_qc/2_cleaned_reads
filtlong --min_length 1000 --keep_percent 90 $WORK/genomics/0_raw_reads/long_reads/*.gz | gzip > $WORK/genomics/2_long_reads_qc/2_cleaned_reads/241155E_cleaned_filtlong.fastq.gz

## 2.3 Nanoplot cleaned
cd $WORK/genomics/2_long_reads_qc/2_cleaned_reads
mkdir -p $WORK/genomics/2_long_reads_qc/3_nanoplot_cleaned
NanoPlot --fastq $WORK/genomics/2_long_reads_qc/2_cleaned_reads/*.gz \
 -o $WORK/genomics/2_long_reads_qc/3_nanoplot_cleaned -t 32 \
 --maxlength 40000 --minlength 1000 --plots kde --format png \
 --N50 --dpi 300 --store --raw --tsv_stats --info_in_report

micromamba deactivate
echo "---------long reads cleaning completed Successfully---------"

module purge
jobinfo
```

# Day 7

## Man muss die nötigen Module laden

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=anvio_pangenomics
#SBATCH --output=anvio_pangenomics.out
#SBATCH --error=anvio_pangenomics.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

# create new folder
mkdir $WORK/pangenomics/02_anvio_pangenomics
```


# Dann muss man die Daten herunterladen

```
curl -L https://ndownloader.figshare.com/files/28965090 -o V_jascida_genomes.tar.gz
tar -zxvf V_jascida_genomes.tar.gz
ls V_jascida_genomes
```

# Dann muss man contigs.dbs aus den .fasta Dateien machen

```
cd $WORK/pangenomics_test/V_jascida_genomes/

ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}_scaffolds.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o V_jascida_${g}.db \
                              --num-threads 4 \
                              -n V_jascida_${g}
done

# annotate contigs.db
for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-scg-taxonomy -c $g --num-threads 4
donecd $WORK/pangenomics_test/V_jascida_genomes/

ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}_scaffolds.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o V_

## Dann können die contigs.db visualisiert werden


### man muss das Terminal öffnen und folgende Befehle eingeben:
```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats $WORK/pangenomics_test/V_jascida_genomes/*db
```
Randbemerkung: Die unteren Schritte wurden nicht ausgeführt


```
srun --reservation=biol217 --pty --mem=16G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8
anvi-display-contigs-stats $WORK/pangenomics_test/V_jascida_genomes/*db
```
### in einem neuen Terminal dann 

```
ssh -L 8060:localhost:8080 sunam238@caucluster.rz.uni-kiel.de

ssh -L 8080:localhost:8080 n100
```

### weitere genomes Dateien keieren

```
cd $WORK/pangenomics_test/V_jascida_genomes
anvi-script-gen-genomes-file --input-dir . -o ./external-genomes.txt
```

## weiterhin muss auf Kontamination geprüft werden

```
cd $WORK/pangenomics_testV_jascida_genomes
anvi-estimate-genome-completeness -e external-genomes.txt
```

# jetzt müssen die Contigs für das Refinement visualisiert werden:

```
cd $WORK/pangenomics_test/V-jascida_genomes/
anvi-profile -c V_jascida_52.db \
             --sample-name V_jascida_52 \
             --output-dir V_jascida_52 \
             --blank
```

create bin V_jascida_52_CLEAN and store it as default

```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-interactive -c V_jascida_52.db \
                 -p V_jascida_52/PROFILE.db
```

Splitting the genome in our good bins

```# Here are the files you created
#V_jascida_52_SPLIT/V_jascida_52_CLEAN/CONTIGS.db
anvi-split -p V_jascida_52/PROFILE.db \
           -c V_jascida_52.db \
           -C default \
           -o V_jascida_52_SPLIT



sed 's/V_jascida_52.db/V_jascida_52_SPLIT\/V_jascida_52_CLEAN\/CONTIGS.db/g' external-genomes.txt > external-genomes-final.txt
```


## Computation of the pangenome

```
anvi-gen-genomes-storage -e external-genomes-final.txt \
                         -o V_jascida-GENOMES.db

anvi-pan-genome -g V_jascida-GENOMES.db \
                --project-name V_jascida \
                --num-threads 4                         
```

## Display of the pangenome

```

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8_biol217

anvi-display-pan -p V_jascida/V_jascida-PAN.db \
                 -g V_jascida-GENOMES.db
```

# Und jetzt muss das ganze für die eigenen ausgesuchten Genome
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=pangenomics_own
#SBATCH --output=pangenmics_own.out
#SBATCH --error=pangenomics_own.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd $WORK/pangenomics/genomes_own_samples/

ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_2.5K.fasta \
                              -o own_samples_${g}.db \
                              --num-threads 4 \
                              -n own_samples_${g}
done

# annotate contigs.db
for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-scg-taxonomy -c $g --num-threads 4
done
```
Dann muss für die Visualisierung im Terminal eingegeben werden:

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats $WORK/pangenomics/genomes_own_samples/*db
```