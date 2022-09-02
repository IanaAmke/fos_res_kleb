#!/bin/bash
#SBATCH --job-name=copyfasta
#SBATCH --account=js66
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1

cp -r /home/iamke/js66/data/klebsiella/pneumoniae/spark_KpSC/assemblies/phenotype /home/iamke/js66_scratch/iana/spark_KpSc_fasta/
