# Metagenomics

Introduction

Metagenomics is the study field of genetic content of microorganisms present in an environment (Offiong et al. 2023). Metagenomics offers insight into the functional gene composition of microbial communities and thus gives a much broader description as previously achieved(Thomas et al. 2012). The cultivation-independent Metagenomics field has enabled profound discoveries regarding the microbial world, offering insights into the diversity and functionalities of microorganisms (Lui et al. 2020; Benz and Mitra, 2023). Its applications span diverse fields including human health, agriculture, and the food industry (Patel et al. 2022). 

Methods

The dataset the Metagenomics analysis bases on is published by Fischer et al. 2018 and includes 16S amplicon sequence analysis. Subsequently, metagenomes were sequenced from the same samples. The samples were obtained from a mesophilic agricultural biogas plant located near Cologne, Germany. Sampling was conducted approximately monthly over a span of 587 days. Fischer et al. reported that the ammonia concentrations during the sampling period exceeded the beneficial ammonia concentration of ~3.0 g/l for acetoclastic methanogenesis, with environmental temperature and pH playing vital roles. 
The objective of this study was to construct Metagenome Assembled Genomes (MAGs) from raw reads through a comprehensive workflow. This process involved preprocessing raw reads, assembling contigs, evaluating assembly quality, binning contigs into MAGs, and assessing MAG completeness, contamination, and strain heterogeneity. Several tools were employed for this purpose, including fastqc, fastp, megahit, samtools, QUAST, Bowtie2, binsanity, MetaBAT2, DASTool, anvi'o, and GUNC (Versions and references are provided in table 1). The analysis was executed on the Linux-based caucluster, leveraging the Slurm batch system for batch computations. Three specific samples from Fischer et al. were selected for analysis. Quality control of raw reads was conducted using fastqc and fastp, followed by genome assembly with megahit. Contig graphs were visualized using Bandage, and assembly quality was assessed using QUAST. Before genome binning, fasta sequence IDs were formatted using anvi'o, and raw reads were mapped onto contigs with Bowtie2. Genome binning was performed using anvi'o with Metabat2 and MaxBin2, followed by quality assessment and estimation of genome completeness. Finally, the resulting data were visualized and evaluated interactively using anvi'o.


Results (Code and figures)

- Quality control of raw reads

In the initial stage of the analysis, the quality of the sequenced data was assessed using fastqc and fastp tools. FastQC provided an overview of basic quality control metrics, particularly focusing on the phred quality score to evaluate the accuracy of base reading. To facilitate organization, a folder was created to store the results using the mkdir command. The analysis began by looping over all files with the .gz extension, representing compressed .fastq files, within the current directory. Alternatively, the second command was utilized for individual files, with both commands directing the output to the designated output_folder.

In our study, fastp was employed to process the reads using various parameters. The implemented command facilitated the processing of paired-end reads by specifying two different inputs for R1 and R2 files within the loop structure, reflecting the complexity inherent in handling paired-end data.

The fastp-processed data was then utilized to execute genome asseblies employing megahit, an NGS assembler known for its efficiency in handling metagenomes coassembly across multiple samples. The assembly was performed with specific parameters tailored for metagenomic data analysis, including a minimum contig length of 1000 base pairs (bp). The resulting contigs were stored in the output folder, with each sample producing an assembly file labeled as final.contigs.fa. It is important to note that these contigs represent only a fraction of the longer context to which they belong.

For visualization of the contig graph in Bandage, the next step involved converting the intermediate contigs in FASTA format to SPAdes-like FASTG format. Subsequently, the generated FASTG file was loaded into Bandage, a GUI-based program, for visualization and analysis. The final.contigs.fastg file served as the input for Bandage, enabling the exploration and interpretation of the contig graph using its graphical user interface (GUI).

![Bandage] (Downloads/Visualization_Contigs_Metagenomics.png)




Discussion

# Folie 16 von Cynthia miteinbeziehen
What can Metagenomics not do?

references for later discussion:
https://doi.org/10.21203/rs.3.rs-1955526/v2 (Impact of microbial genome completeness on functional metagenomics)



# references 
Benz, S., & Mitra, S. (2023). From Genomics to Metagenomics in the Era of Recent Sequencing Technologies. In Metagenomic Data Analysis (pp. 1-20). New York, NY: Springer US.

Fischer, M. A., Güllert, S., Refai, S., Künzel, S., Deppenmeier, U., Streit, W. R., & Schmitz, R. A. (2019). Long‐term investigation of microbial community composition and transcription patterns in a biogas plant undergoing ammonia crisis. Microbial Biotechnology, 12(2), 305-323.

Lui, L. M., Nielsen, T. N., & Arkin, A. P. (2020). A method for achieving complete microbial genomes and better quality bins from metagenomics data. bioRxiv.

Offiong, N. A. O., Edet, J. B., Shaibu, S. E., Akan, N. E., Atakpa, E. O., Sanganyado, E., ... & Okoh, A. (2023). Metagenomics: an emerging tool for the chemistry of environmental remediation. Frontiers in Environmental Chemistry, 4, 7.

Patel, T., Chaudhari, H. G., Prajapati, V., Patel, S., Mehta, V., & Soni, N. (2022). A brief account on enzyme mining using metagenomic approach. Frontiers in Systems Biology, 2, 1046230.

Thomas, T., Gilbert, J., & Meyer, F. (2012). Metagenomics-a guide from sampling to data analysis. Microbial informatics and experimentation, 2, 1-12.

