import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header = re.findall(r'>[\w\.-]+[.][\w\.-]+Fos[A-Z,0-9]|^[>Fos|>sp|>pdb][\w\.-]+|^[>sp|>pdb][\w\.-]+Fos[A-Z,0-9]', line)[0]
            seq_dict[header] = ""
        else:
            seq_dict[header] = seq_dict[header] + line
        print(seq_dict[header])
            
            
