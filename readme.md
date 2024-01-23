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

# Task: How to add links and images to the end file


##  code for Sequencing files
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
