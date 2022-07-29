import re

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header = re.findall(r'>[\w\.-]+[.][\w\.-]+Fos[A-Z,0-9]', line)
            header_2 = re.findall(r'>Fos|sp|pdb', line)
            #print(header)
            print(header_2)
            
