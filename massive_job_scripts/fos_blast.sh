#!/bin/bash
#SBATCH --job-name=fos_blast
#SBATCH --account=js66
#SBATCH --time=01:00:00
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1

# load modules on M3
module load blast/2.2.30
module load python/2.7.15-gcc5

# blast Kp genomes from samples against fos genes
python /home/iamke/js66/holtlab/geneBlastCheck.py -s /home/iamke/js66_scratch/iana/fos_genes_all_db_clean.fasta /home/iamke/js66_scratch/iana/Egli_MIC_fasta/*.fasta > Egli_MIC_gene_hits.txt

python /home/iamke/js66/holtlab/geneBlastCheck.py -s /home/iamke/js66_scratch/iana/fos_genes_all_db_clean.fasta /home/iamke/js66_scratch/iana/Huynh2020_fasta/*.fasta > Huynh2020_gene_hits.txt

python /home/iamke/js66/holtlab/geneBlastCheck.py -s /home/iamke/js66_scratch/iana/fos_genes_all_db_clean.fasta /home/iamke/js66_scratch/iana/sands2021_arthropods_fasta/*.fasta > sands2021_arthropods_gene_hits.txt

python /home/iamke/js66/holtlab/geneBlastCheck.py -s /home/iamke/js66_scratch/iana/fos_genes_all_db_clean.fasta /home/iamke/js66_scratch/iana/sands2021_BARNARDS_fasta/*.fasta > sands2021_BARNARDS_gene_hits.txt

python /home/iamke/js66/holtlab/geneBlastCheck.py -s /home/iamke/js66_scratch/iana/fos_genes_all_db_clean.fasta /home/iamke/js66_scratch/iana/spark_KpSc_fasta/*.fasta > spark_KpSc_gene_hits.txt

