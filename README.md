# Fosfomycin Resistance in *K. pneumoniae* 
## Project summary
*K. pneumoniae* is a Gram-negative bacterium that causes opportunistic and community-acquired infections. *K. pneumoniae* easily develops resistance to multiple antibiotics as it can acquire numerous antibiotic resistance genes (ARGs) from the environment and from other bacteria. As a result, there exists multidrug resistant (MDR) and extremely resistant (XDR) strains which are highly resistant to almost all the antibiotics currently in use.

Fosfomycin is an antibiotic that has been in use in the treatment of uncomplicated UTIs and is currently being trialled, either as a monotherapy or in combination with other drugs, in the treatment of MDR. However, mechanisms of resistance have been reported with one of them being the presence of *fos* genes that encode for Fos enzymes which inactivate fosfomycin through the nucleophilic addition of substrates such as glutathione, bacillithiol, and water to the first carbon of its epoxide ring.  

Mechanisms of fosfomycin resistance have been well-studied in bacteria such as *E. coli* where they have been shown to contribute to higher minimum inhibitory concentration (MIC) measurements and worse clinical outcomes for treatment with fosfomycin. However, there have been a few studies on mechanisms of resistance in *K. pneumoniae* and only a few have gone into depth on the impact of the presence of mutations in the fosfomycin transport systems and drug target and genes encoding fosfomycin-modifying enzymes.

Therefore, the main aim of this project was to study fosfomycin resistance mechanisms in *K. pneumoniae* and in particular, the impact of the presence of *fos* genes to fosfomycin resistance levels in the bacteria.  

## Data analysis
The data used for this project was made up of *K. pneumoniae* assembled genomes that were matched with antibiogram data. *fos* genes used (data/final_seq_data/fos_genes_final_clean.fasta) were curated from the following public databases: CARD, ResFinder, Kleborate, and GenBank, and were queried against the sample genomes to find which *K. pneumoniae* genomes contained *fos* gene sequences. *fos* genes obtained from the different databases are contained in the folder seq_data (data/seq_data/). 

Some of the tools used for analysis include:
* R for data visualisation and phylogenetic tree annotation
  + ggplot2 (v3.3.6) to create bar and violin plots
  + ggupset (v0.3.0) to create UpSet plots
  + ggtree (v3.4.2) for phylogenetic tree annotation
  + ggtreeExtra (v1.6.0) to create heatmaps representing tree metadata
* Command-line tools for data cleaning and phylogenetic analysis
  + CD-HIT (v4.8.1) to delete duplicate sequences
  + SeqKit (v2.2.0) to delete duplicates and remove gaps from FASTA files
  + BLAST (v2.2.30) to query *K. pneumoniae* sample genomes against *fos* genes
  + MAFFT (v7.505) to align *fos* gene sequences extracted from the sample genomes
  + FastTree (v2.1.11) for phylogenetic analysis - maximum likelihood phylogenetic trees
* Python for FASTA file quality control and analysis. Python scripts, written by members of the Holt lab, were used to automate some of the analysis on the Massive M3  high performance computing system (HPC).

## R and Python scripts written for this project
### R scripts
* summ_antibiogram.R was used to create bar plot summaries of the datasets being used for analysis
* fos_blast_uiquehits.R was used to create bar plots visualising the distribution of genomes with *fos* genes across measurements of the different susceptibility testing methods i.e. MIC (broth dilution) and disk diffusion. It also gives an output of violin + UpSet plots showing the distribution of genomes with *fos* genes across measurements of the different susceptibility testing methods in addition to the number of genomes with specific *fos* genes and *fos* gene combinations for genomes with more than one *fos* gene
* fos_gene_tree.R and fos_genome_tree.R were used to create a gene and genome level tree respectively with heatmaps showing gene copy number and different AST measurements 
* fos_genome_tree_circular.R was used to create a circular genome tree with a heatmap showing gene copy number
* mutate_res_phenotype_func.R was a function created to automatically fill the phenotype of a sample based on EUCAST/CLSI AST breakpoints 

### Python scripts
* fasta_header_edit_AA.py and fasta_header_edit_gene.py were used to edit the FASTA file header to only include the gene name and ID
* filter_fosgenes.py was used to check whether all *fos* genes had been extracted from the blast hits (the most common *fos* gene was used as a query sequence used to extract *fos* gene seqeunces from the sample genomes) 
* start_methionine_AA.py was used to check whether Fos amino acid sequences started with a methionine
