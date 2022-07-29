import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'\_([A-Z,0-9,]*)', line)
            #header2 = re.findall(r'>[\w\.-]+[.][\w\.-]+Fos[A-Z,0-9]|^[>Fos|>sp|>pdb][\w\.-]+|^[>sp|>pdb][\w\.-]+Fos[A-Z,0-9]', line)[0]
            print(header1)
            