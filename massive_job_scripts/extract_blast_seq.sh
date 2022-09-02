#!/bin/bash
#SBATCH --job-name=extract_seq
#SBATCH --account=js66
#SBATCH --time=01:00:00
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1

# load modules on M3
module load blast/2.2.30
module load python/2.7.15-gcc5

# extract alignment of BLAST hits from a set of input seq
python /projects/js66/holtlab/kats_script_folder/blast2Align.py -s query_seq.fasta -l 80 -i 80 -o Egli_MIC_geneHits.mfasta /home/iamke/js66_scratch/iana/Egli_MIC_fasta/*.fasta
python /projects/js66/holtlab/kats_script_folder/blast2Align.py -s query_seq.fasta -l 80 -i 80 -o Huynh2020_geneHits.mfasta /home/iamke/js66_scratch/iana/Huynh2020_fasta/*.fasta
python /projects/js66/holtlab/kats_script_folder/blast2Align.py -s query_seq.fasta -l 80 -i 80 -o sands2021_arthropods_geneHits.mfasta /home/iamke/js66_scratch/iana/sands2021_arthropods_fasta/*.fasta
python /projects/js66/holtlab/kats_script_folder/blast2Align.py -s query_seq.fasta -l 80 -i 80 -o sands2021_BARNARDS_geneHits.mfasta /home/iamke/js66_scratch/iana/sands2021_BARNARDS_fasta/*.fasta
python /projects/js66/holtlab/kats_script_folder/blast2Align.py -s query_seq.fasta -l 80 -i 80 -o spark_KpSc_geneHits.mfasta /home/iamke/js66_scratch/iana/spark_KpSc_fasta/*.fasta
