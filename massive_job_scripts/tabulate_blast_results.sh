#!/bin/bash
#SBATCH --job-name=blast_tabulate
#SBATCH --account=js66
#SBATCH --time=01:00:00
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1

# load modules
module load python/2.7.15-gcc5

# tabulate raw blast output
python /projects/js66/holtlab/tabulateBlastHits.py -i Egli_MIC_gene_hits.txt -o Egli_MIC_gene_search_results
python /projects/js66/holtlab/tabulateBlastHits.py -i Huynh2020_gene_hits.txt -o Huynh2020_gene_search_results
python /projects/js66/holtlab/tabulateBlastHits.py -i sands2021_arthropods_gene_hits.txt -o sands2021_arthropods_gene_search_results
python /projects/js66/holtlab/tabulateBlastHits.py -i sands2021_BARNARDS_gene_hits.txt -o sands2021_BARNARDS_gene_search_results
python /projects/js66/holtlab/tabulateBlastHits.py -i spark_KpSc_gene_hits.txt -o spark_KpSc_gene_search_results
