import re

seq_dict = {}

with open("fos_AA_final_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header = re.split('\s', line)[1]
        