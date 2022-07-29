import re

seq_dict = {}

with open("fos_genes_final.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'fos[A-Z,0-9]*', line)[0]
            header2 = re.findall(r'\_([\w\.-]+)', line)
            print(header2)