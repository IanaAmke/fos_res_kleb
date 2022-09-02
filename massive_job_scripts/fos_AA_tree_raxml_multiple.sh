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
raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned.phy -n fos_AA_ML_1 -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 12345

raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned.phy -n fos_AA_ML_2 -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 12345

raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned.phy -n fos_AA_ML_3 -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 12345

raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned.phy -n fos_AA_ML_4 -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 12345

raxmlHPC-PTHREADS-SSE3 -T 8 -s /home/iamke/js66_scratch/iana/fos_AA_with_outliers_aligned.phy -n fos_AA_ML-5 -f a -m PROTGAMMAAUTO -x 12345 -N 1000 -p 12345

