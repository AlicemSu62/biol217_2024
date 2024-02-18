# Pangenomics

## Introduction

The pangenome serves as a tool for exploratory analyses of genome functions and functional relationships. The field of pangenomics compares complete genomes or HQ MAGs. Hereby, the functional enrichment based on environmental conditions and evolutions get calculated and the average nucleotide idendity (ANI) scores between genomes can be computed and visualized.

Pangenomics used to compare the genetic relatedness including the core genome, the accessory genes as well as the singletons. Today, Pangenomics compares closely related HQ genomes and clusters genes based on the amino acid (AA) sequence comparison. General Quality information given in pangenomics analysis are singletons, redundancy, completion and, the GC-content. THe initial gene-cluster analysis starts with calculating the similarity of each AA-sequence between each genome, which is the Diamond Blast method (Step 1). Then, weak similarity hits such as incomplete and trunctated genes are being removed (Step 2). After removing those hits, clustering into groups based on Marcov Cluster Algorithm (MCL) follows (Step 3). Linking to the algorithm is a depictins as a dendrogram and rings to organize the pangenome (Step 4). In the Gene-frequency analysis, 3 main genome types are distinguished. The core genome are genes shared between all genomes and are often essential metabolic genes. One gene in the core genome are the single copy core genes. The accessory genome are genes shared by some genomes and they are aquired through environmental adaptation. Unique genes, often referred to as singletons, are genes present only in one genome and mobile genetic elements such as Phages, transposons and CRISPRs.

Single copy core genes (SCG)are essential genes encoded once per genome and shared among all genomes in a dataset. SCGs are utilized for the assessment of MAG purity and contamination, also for the genome completeness and phylogenomics/phylogenetic trees. 

Geometric and funtional homogeneity are used to detect highly conserved or variable regions. Homogeneity is an indicator for evolutionary changes within a genomic region. Geometric homogeneity is based on the configuration of genes within a gene cluster, considers gaps and residues structural configurations. Functional homogeneity on the other hand is based on biochemical properties across genes clusters and considers changes in AA residues which affect the biochemistry of a protein. 

In order to figure out geometric homogenity, AA sequences are aligned (Step 1), Gap/residue patterns are determined (Step 2), and site- and gene-level changes of residues combined (Step3). (All in all, the order of gene cluster plus the residue changes equals geometric homogeneity)

In funtional homogeneity, AAs are categorized into biochemical groups based on their side chains. Therefore, the algorithm is considering, if each residue is conserved or changed and if changed, it will ask if the new residue will likely change conformation based on its polarity/size.

Functional annotation is based on known information of genes/proteins and used to decorate genomes with information and on AA sequences. Thereby, it accessess available databases such as KEGG/KOfams, COG and Pfam (explained in a bit). Annotations are using Hidden Markov Models (HMM) are statistical models based on a scoring system and detect homology based on AA sequence. HMMs can model protein sequences depending on the features of the protein. In a reference, HMMs can be described as types of directed graphs. Databases mentioned previously are the Kyoto Encyclopedia of Genes and Genomes (KEGG Japan) and the Clusters of Orthologous Genes (COG NCBI). Both are libraries of genes functions (characterized orthologues). With the help o the data bases whole pathways can be visualized and metabolic modules can be assigned. Both databases are based on experimental results and available as online resources. 

the InterPro database provides an integrative classification of protein sequences into families, and identifies functionally important domains and conserved sites. It also provides the prediction of proteins and domains and combines classifications and approaches from member databases (HMMs, patterns, fingerprints)

AlphaFold is an Ai-based system developed by DeepMind. It predicts 3D Protein structures from exported AA sequences and enables functional interpretation and evolution experiments and gives information complementary to pangenome analysis.

For pangenomics, it is important to know when to cluster what. Gene clusters are grounded on the AA similarities between genomes. After clustering the genes, gene frequencies are being analyzed and with the allocation of genes into the core genome, accessory genome and Singletons questions as "which genes are essential" , and "which genes are individual" can be answered. Geometric and functional homogeneity look for changes in gene cluster patterns or biochemical function and give information on evolution analysis through region with high/low variability. at last, with the functional annotations, it is looked for specific metabolic pathways and features and asked, what the organism can do.

Pangenomics results



## Pangenomics - comparing genomes using ANVIO

(the following must be done on the caucluster, because I have not downloaded the required modules yet)



