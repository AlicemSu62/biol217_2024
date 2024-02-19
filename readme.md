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

the caucluster can be accessed via command `ssh -X sunam238@caucluster.rz.uni-kiel.de`




### Task: How to add links and images to the end file

```
![Picture](path/to/picture)

```
-------

# Day 2

##  code for Sequencing files

zuerst muss gesagt werden, was der Umstand für das Programmieren ist. Wir haben in dem Fall jetzt sequence files, 6 Stück insgesamt, 3 davon sind normal und 3 davon reversed reads. Am Ende ergeben 1 normal read und ein reversed read paired end readfiles. Die Files habe den Datentypen fastq.gz.
Mit fastqc (qc für QualityControl) verwirklicht man eine Visualisierung der fastq Dateien. HTMLs werden ausgegeben, die man dann in einen Internetbrowser kopieren kann und die visualisierte Datei sieht. 

In the initial stage of the analysis, the quality of the sequenced data (6 reads, 3 forward and 3 reverse - 1 of forward and reverse make one paired end readfile - so in sum, we get 3 paired end readfiles) was assessed using fastqc and fastp tools. Fastqc provided an overview of basic quality control metrics, particularly focusing on the phred quality score to evaluate the accuracy of base reading. To facilitate organization, a folder was created to store the results using the mkdir command. The analysis began by looping over all files with the .gz extension, representing compressed .fastq files, within the current directory. Alternatively, the second command was utilized for individual files, with both commands directing the output to the designated output_folder.

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

jobs can be submitted through the command sbatch

```
sbatch jobscript.sh
```


the status of submission can be checked via squeue (if job is submitted, if job is still processing)

```
squeue -u sunam238
```
Die html s auf den eigenen PC kopieren, dafür einfach neues Terminal-Fenster öffnen, denn das aktuelle Terminal greift auf das cau-cluster zu.

```
scp sunam238@caucluster.rz.uni-kiel.de:/work_beegfs/suname238/Metagenomics/1_fastqc/*.html.
```

## fast p

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

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8


fastp -i BGR_130305_mapped_R1.fastq.gz  -I BGR_130305_mapped_R2.fastq.gz  -R BGR_130305_reportfile -o ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130305_mapped_R2_clean.fastq.gz -t 6 -q 20

fastp -i BGR_130527_mapped_R1.fastq.gz  -I BGR_130527_mapped_R2.fastq.gz  -R BGR_130527_reportfile -o ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -t 6 -q 20

fastp -i BGR_130708_mapped_R1.fastq.gz  -I BGR_130708_mapped_R2.fastq.gz  -R BGR_130708_reportfile -o ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -O ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -t 6 -q 20

```
---------------------

## megahit

The fastp-processed data was then utilized to execute genome asseblies employing megahit, an NGS assembler known for its efficiency in handling metagenomes coassembly across multiple samples. The assembly was performed with specific parameters tailored for metagenomic data analysis, including a minimum contig length of 1000 base pairs (bp). The resulting contigs were stored in the output folder, with each sample producing an assembly file labeled as final.contigs.fa. It is important to note that these contigs represent only a fraction of the longer context to which they belong.

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

For visualization of the contig graph in Bandage, the next step involved converting the intermediate contigs in FASTA format to SPAdes-like FASTG format. Subsequently, the generated FASTG file was loaded into Bandage, a GUI-based program, for visualization and analysis. The final.contigs.fastg file served as the input for Bandage, enabling the exploration and interpretation of the contig graph using its graphical user interface (GUI).

```
megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg                   
```

Questions
Please submit your generated figure and explain in your own words what you can see (keep it short).

***The graph shows all contigs resulting from megahit with a total of 57414 nodes, 8 edges and a total base length of 145675865 kbp. Other data given in Bandage are N50=2963 bp.***



---------------------------

# Day 3

Quast, a QUality ASsessment Tool, is utilized to evaluate genome assembly, particularly the results obtained from megahit. To initiate the evaluation process, navigate to the Quast folder using the terminal and execute the following jobscript:

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

cd /work_beegfs/sunam238/Metagenomics

metaquast -t 6 -o ./3_metaquast -m 1000 ./3_coassembly/final.contigs.fa
```

----

html-file can be downloaded or following command can be given into the local terminal:
```
scp sunam238@caucluster.rz.uni-kiel.de:/work_beegfs/sunam238/Metagenomics/3_metaquast/report.html 
```

Questions
What is your N50 value? Why is this value relevant?
How many contigs are assembled?
What is the total length of the contigs?


***N50 value and relevance: N50 where the lengths of aligned blocks are counted instead of the contig lengths. I.e., if a contig has a missassembly with respect to the reference, the contig is broken into smaller pieces. This metric is computed only if a reference genome is computed. ***
How many contigs are assembled: 57414 
- Methanoculleus_bourgensis_MS2:	958
- Porphyromonadaceae_bacterium_ING2_E5B:	399
- not_aligned:	56057
What is the total length of the contigs: 145675865 kbp


## Binning 

For Binnning, the initial step involved formatting fasta sequence IDs using anvi'o to ensure subsequent functionality during sequence matching and mapping processes.

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

cd /work_beegfs/sunam238/Metagenomics/

anvi-script-reformat-fasta ./3_coassembly/final.contigs.fa -o /work_beegfs/sunam238/Metagenomics/3_coassembly/contigs.anvio.fa --min-len 1000 --simplify-names --report-file name_conversion.txt
```

Following this, clean reads were mapped onto assembled contigs using Bowtie2, with prior indexing of the mapping reference fasta file to expedite the mapping process.

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


cd work_beegfs/sunam238/Metagenomics/3_coassembly/

module load bowtie2
bowtie2-build contigs.anvio.fa contigs.anvio.fa.index

bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130305_mapped_R2_clean.fastq.gz -S map130305.sam

bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -S map130527.sam

bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -S map130708.sam
```  

The SAM files are converted into a binary alignment and map (BAM) file with the extension .bam using samtools.

```
module load samtools
samtools view -bS sam_file.sam > bam_file.bam
```
or in a loop with the adequate path: `for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done`

----------

## Mapping

Anvi'o was then employed to preprocess contigs data, involving computations of k-mer frequencies, soft-splitting of contigs longer than 20,000 bp, and identification of open reading frames using Prodigal.
Subsequent Hidden Markov Model (HMM) searches on contigs were performed using anvi'o to identify genes with known functions, leveraging multiple default bacterial single-copy core gene collections to aid in gene hit identification.
Upon preparation of the contigs database and optional execution of HMM searches, anvi-display-contigs-stats facilitated a rapid assessment of assembly output and an estimation of recoverable bacterial and archaeal genomes.
Genome binning was executed using Anvi'o, encompassing steps such as sorting and indexing of .bam files and establishment of an Anvi'o profile to store sample-specific contig information. 


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=anvi
#SBATCH --output=anvi.out
#SBATCH --error=anvi.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam238/Metagenomics/3_coassembly/


anvi-gen-contigs-database -f contigs.anvio.fa -o contigs.db -n 'biol217'


anvi-run-hmms -c contigs.db


srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --nodelist=node002 /bin/bash
```

An anvi'o profile retains sample-specific details regarding contigs. By profiling a BAM file with anvi'o using the command anvi-profile, a unified profile is generated. Anvi-profile further processed each contig, providing insights into mean coverage, standard deviation of coverage, and single-nucleotide variants.

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

anvi-profile -i map130708.sam..bam -c contigs.db --output-dir OUTPUT_DIR
mkdir /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/profiling/
for i in `ls *.sorted.bam | cut -d "." -f 1`; do anvi-profile -i "$i".bam.sorted.bam -c work_beegfs/sunam238/Metagenomics/3_coassembly/contigs.db -o /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/profiling/”$i”; done
```

important note: Cynthia sent me required data because of previous errors which are in the folder CC `/work_beegfs/sunam238/Metagenomics/CC/`
----



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

anvi-interactive -p /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam238/Metagenomics/CC/4_mapping/contigs.db -C YOUR_COLLECTION
```
-----

# Day 4 (but really Day 3 continued)


Day 3 continued


anvi-display provides basic statistics of the contigs database, aiding in the assessment of the assembly output and enabling estimation of the number of bacterial and archaeal genomes that can be recovered.

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

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


# Binning with METABAT and MAXBIN2

Finally, genome binning was performed utilizing Metabat2 and MaxBin2.



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
Number of Archaea bins from MetaBAT2: 3
Number of Archaea bins from Maxbin2: 2


### MAGs Quality Estimation

The completeness and contamination levels of the genomes were estimated to assess the quality of the bins. This evaluation was performed using anvi-estimate-genome-completeness.

```
cd /work_beegfs/sunam238/Metagenomics/
anvi-estimate-genome-completeness -c ./CC/4_mapping/contigs.db -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -C METABAT2
```

```
cd /work_beegfs/sunam238/Metagenomics/
anvi-estimate-genome-completeness -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db --list-collections
```


The data required for this assessment was already available in the HTML file generated from binning.

Subsequently, the results were visualized and evaluated using anvi-interactive. This tool allowed for manual inspection and manipulation of bins.
```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-interactive -p ./CC/5_anvio_profiles/merged_profiles/PROFILE.db -c ./CC/4_mapping/contigs.db -C YOUR_COLLECTION
```
Which binning strategy gives the best quality for the ARCHAEA bins?
***METABAT works with 375 contigs creating 48 bins while MAXBIN works with 1381 contigs and makes 69 bins. Out of these 69 bins, one is an Archaea bin. METABAT on the other hand, has 3 Archaea bins. In order to say, which bin gives better quality Archaea bins, the parameter Completion has to be looked over for the Archaea bins for each binning type. MAXBIN_10 has a completion of 96.05%, METABAT_14, METABAT_8 and METABAT_39 have completions of 43.68%, 97.37%, and 50.00%.  ***
How many Archaea bins does one get that are of high quality? How many Bacteria bins does one get that are of High Quality?
***Bins that are of high quality are those with a completion over 90% percent and a contamination less than 5%. METABAT_8 would be the only high quality bin with 97.37% completion and about 5% redundancy.***

Archaea bins found:

MAXBIN_10 : Compl. 96.05%} from displayed summary
MAXBIN_10: Compl. 96.05%
MAXBIN_22: Compl. 0%



METABAT_8 : Compl. 97.37% 
METABAT_14 : Compl. 43.68%
METABAT_39: Compl. 50 %} from displayed summary
METABAT_19 : 97.37%
METABAT_23 : 39.47%
METABAT_41 : 50.00%

# Day 4 ()

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
cd /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/ARCHAEA_BIN_REFINEMENT

mkdir GUNC

for i in *.fa; do gunc run -i "$i" -r /work_beegfs/sunam238/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC/"i" --threads 10 --detailed_output; done


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
anvi-refine -c /work_beegfs/sunam238/Metagenomics/CC/4_mapping/contigs.db -C METABAT -p /work_beegfs/sunam238/Metagenomics/CC/5_anvio_profiles/merged_profiles/PROFILE.db --bin-id METABAT__44
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

anvi-script-get-coverage-from-bam -b /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/map130305.sam.bam.sorted.bam -c /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/profiling/map130305/AUXILIARY-DATA.db/ -C YOUR_COLLECTION -m contig -o /work_beegfs/sunam238/Metagenomics/CC/

anvi-script-get-coverage-from-bam -b /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/map130527.sam.bam.sorted.bam 

anvi-script-get-coverage-from-bam -b /work_beegfs/sunam238/Metagenomics/5_anvio_profiles/map130708.sam.bam.sorted.bam 

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

ls *fasta > genomes.txt

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

# Create external genomes file

```
anvi-script-gen-genomes-file --input-dir $WORK/pangenomics/genomes_own_samples -o external-genomes.txt
```

Auf Kontamination prüfen

```
cd $WORK/pangenomics/genomes_own_samples
anvi-estimate-genome-completeness -e external-genomes.txt
```


Jetzt muss das pangenome computisiert werden

```
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o $WORK/pangenomics/genomes_own_samples/own_sample-GENOMES.db

anvi-pan-genome -g own_sample-GENOMES.db \
                --project-name ownsample \
                --num-threads 4 --enforce-hierarchical-clustering
```

zuletzt muss das Pangenom visualisiert werden

```


module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8_biol217

anvi-display-pan -p $WORK/pangenomics/genomes_own_samples/own_sample-GENOMES.db \
                 -g own_sample-GENOMES.db
```

funktioniert leider nicht. Was ich jetzt versuche, ist, dass ich erstmal nur 5 Genome visualisiere, da die Anzahl der Gene Grund für den Error sein könnte:

ich muss nochmal von vorne anfangen und die Genome in zwei separate Ordner splitten:

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

cd $WORK/pangenomics/pangenome_new/

# ls *fasta > genomes.txt

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

# Day8

zuerst müssen Module miniconda und READemption geladen und ein Project_Path erstellt werden:

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

export http_proxy=http://relay:3128
export https_proxy=http://relay:3128
export ftp_proxy=http://relay:3128

cd $WORK/RNAseq/
conda activate reademption
reademption create --project_path READemption_analysis --species salmonella="Salmonella Typhimurium"
```

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-04:00:00
#SBATCH --job-name=reademption_tutorial
#SBATCH --output=reademption_tutorial.out
#SBATCH --error=reademption_tutorial.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna

sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa
wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz

gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2

reademption align -p 4 --poly_a_clipping --project_path READemption_analysis

reademption coverage -p 4 --project_path READemption_analysis

reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 -c InSPI2,InSPI2,LSP,LSP -r 1,2,1,2 --libs_by_species salmonella=InSPI2_R1,InSPI2_R2,LSP_R1,LSP_R2 --project_path READemption_analysis

reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
conda deactivate
module purge
jobinfo
```

Open the paper from this Prasse et al. 2017, find out the SRR numbers, quantity of samples and treatments, and write down here:


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-04:00:00
#SBATCH --job-name=reademption_tutorial
#SBATCH --output=reademption_tutorial.out
#SBATCH --error=reademption_tutorial.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load micromamba/1.4.2
micromamba activate 10_grabseqs

export http_proxy=http://relay:3128
export https_proxy=http://relay:3128
export ftp_proxy=http://relay:3128


mkdir /$WORK/RNAseq/READemption_analysis_my/fastq
cd /$WORK/RNAseq/READemption_analysis_my/fastq/

grabseqs -t 4 -m ./metadata.csv SRR4018517
grabseqs -t 4 -m ./metadata.csv SRR4018516
grabseqs -t 4 -m ./metadata.csv SRR4018515
grabseqs -t 4 -m ./metadata.csv SRR4018514

micromamba deactivate

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate reademption

reademption create --project_path READemption_analysis_my --species methanosarcina="Methanosarcina mazei Go1"
```

```
mkdir ../qc_reports
for i in *.fastq.gz; do fastqc -t 4 -o ../qc_reports/fastqc_output $i; done
```

Quality Control

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-04:00:00
#SBATCH --job-name=qc_reademption
#SBATCH --output=qc_reademption.out
#SBATCH --error=qc_reademption.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load micromamba/1.4.2
micromamba activate 01_short_reads_qc

mkdir $WORK/RNAseq/fastq/qc_reports
for i in *.fastq.gz; do fastqc -t 4 -o $WORK/RNAseq/fastq/qc_reports/fastqc_output $i; done
```

sed -i "s/>/>NC_003901.1 /" $WORK/RNAseq/READemption_analysis_my/input/methanosarcina_reference_sequences/methanosarcina_mazei_go1.fna

```
#!/bin/bash
#SBATCH --job-name=reademption_
#SBATCH --output=reademption.out
#SBATCH --error=reademption.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --partition=base
#SBATCH --export=NONE
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate reademption

reademption align -p 4 --poly_a_clipping --project_path READemption_analysis
reademption coverage -p 4 --project_path READemption_analysis
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l wt_1.fasta.gz,wt_2.fasta.gz,mut_1.fasta.gz,mut_2.fasta.gz -c wt,wt,mut,mut -r 1,2,1,2 --libs_by_species methanosarcina=wt_1,wt_2,mut_1,mut_2 --project_path READemption_analysis
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
conda deactivate
jobinfo
```

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-04:00:00
#SBATCH --job-name=rna_seq_methanosarcina
#SBATCH --output=rna_seq_methanosarcina.out
#SBATCH --error=rna_seq_methanosarcina.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate reademption

#2- copy the sequences and files in respective directories
# download the sequences from the NCBI database or github folder named "genome_input"

#3- Processing and aligning the reads
reademption align --project_path READemption_analysis \
	--processes 32 --segemehl_accuracy 95 \
	--poly_a_clipping \
	--fastq --min_phred_score 25 \
	--progress

#4- Coverage
reademption coverage --project_path READemption_analysis \
	--processes 32

#5- Performing gene wise quantification
reademption gene_quanti --project_path READemption_analysis \
	--processes 32 --features CDS,tRNA,rRNA 

#6- Performing differential gene expression analysis 

####NOTE:: Change the names according to your file names in the READemption_analysis/input/reads/ directory
reademption deseq --project_path READemption_analysis \
	--libs mut_1.fastq.gz,mut_2.fastq.gz,wt_1.fastq.gz,wt_2.fastq.gz \
	--conditions mut,mut,wt,wt --replicates 1,2,1,2 \
	--libs_by_species metanosarcina=mut_1,mut_2,wt_1,wt_2

#7- Create plots 
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis

#The whole command will take around 2 hours to run.
conda deactivate
module purge
jobinfo
```

what is tnoar?
TNOAR is the total number of aligned reads. 


aLSO I HAVE TO WRITE DOWN THE EQUATIONS FOR RKPM ET CETERA


Now we have to write down 5 upregulated and 5 downregulated genes in the deseq annotation file to work with them in R

Upregulated
RS_02180
RS_02465
RS_02525
RS_02650
RS_02950

Downregulated
RS_13050
RS_12835
RS_19460
RS_17565
RS_16165

#How to create an image:

![igb](Path/to/image)