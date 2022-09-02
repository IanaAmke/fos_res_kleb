#!/bin/bash
#SBATCH --job-name=sort_contigs
#SBATCH --account=js66
#SBATCH --time=01:00:00
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1

# load modules on M3
module load python/2.7.15-gcc5

# assign specific alleles to contig regions
python /projects/js66/holtlab/sortContigsFromResFinderBlast.py Egli_MIC_gene_hits.txt -t restable > Egli_MIC_uniqueHits.txt
python /projects/js66/holtlab/sortContigsFromResFinderBlast.py Huynh2020_gene_hits.txt -t restable > Huynh2020_uniqueHits.txt
python /projects/js66/holtlab/sortContigsFromResFinderBlast.py sands2021_arthropods_gene_hits.txt -t restable > sands2021_arthropods_uniqueHits.txt
python /projects/js66/holtlab/sortContigsFromResFinderBlast.py sands2021_BARNARDS_gene_hits.txt -t restable > sands2021_BARNARDS_uniqueHits.txt
python /projects/js66/holtlab/sortContigsFromResFinderBlast.py spark_KpSc_gene_hits.txt -t restable > spark_KpSc_uniqueHits.txt
