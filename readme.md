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

```
cd /PATH/TO/SUMMARY/bin_by_bin

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

for i in *.fa; do gunc run -i ? -r /home/sunam226/Databases/gunc_db_progenomes2.1.dmnd --out_dir ? --threads 10 --detailed_output; done

```

```

```