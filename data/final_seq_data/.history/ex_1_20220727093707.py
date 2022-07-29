import re

seq_dict = {}

with open("fos_AA_with_outliers_aligned_trimmed.fasta.reduced") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        header = re.split('\s', line)[1]
        print(header)
        