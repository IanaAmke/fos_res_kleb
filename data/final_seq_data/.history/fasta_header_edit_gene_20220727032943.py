import re

seq_dict = {}

with open("fos_genes_final.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'Fos[A-Z,0-9]*', line)
            print(header1)