import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'Fos[A-Z,0-9]*|sp[A-Z,0-9]|pdb[A-Z,0-9]|Orf[A-Z,0-9]', line)
            header2 = re.findall(r'\_([\w\.-]*)', line)
            print(header1)