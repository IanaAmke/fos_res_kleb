import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'\_([\w\.-]*)', line)
            header2 = re.findall(r'\_(Fos[A-Z,0-9]*)', line)
            print(header1)