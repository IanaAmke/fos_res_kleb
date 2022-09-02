#!/bin/bash
#SBATCH --job-name=fos_ML
#SBATCH --account=js66
#SBATCH --partition=m3h,m3g
#SBATCH --qos=normal
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=24576

# Job modules
module purge
module load raxml/8.2.9

# run RAxML
raxmlHPC-PTHREADS-SSE3 -T 10 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned_trimmed.fasta.reduced -n fos_AA_ML -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 1234

