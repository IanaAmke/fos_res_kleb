#!/bin/bash
#SBATCH --job-name=fos_ML
#SBATCH --partition=m3i,m3m,comp
#SBATCH --account=js66
#SBATCH --qos=normal
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=24576

# Job modules
module purge
module load raxml/8.2.9

# run RAxML
raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_genes_final_clean.fasta -n fos_AA_ML -f a -m GTRGAMMA -x 12345 -N 1000 -p 12345

