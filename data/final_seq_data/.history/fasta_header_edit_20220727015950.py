import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'Fos[A-Z,0-9]*|Orf[A-Z,0-9]', line)[0]
            header2 = re.split('\s', line)[0]
            if line.startswith(">Fos"):
                header2 = re.split('\s', line)[1]
            elif line.startswith(">Orf"):
                header2 = re.split('\s', line)[1]
            print(header1)