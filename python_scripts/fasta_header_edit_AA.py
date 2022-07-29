import re

seq_dict = {}

with open("fos_AA_nogaps_ex.fasta") as fosAA_file:
    for line in fosAA_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'Fos[A-Z,0-9]*|Orf[A-Z,0-9]', line)[0]
            header2 = re.split('\s', line)[1]
            if line.startswith("> Fos"):
                header2 = re.split('\s', line)[2]
            elif line.startswith("> Orf"):
                header2 = re.split('\s', line)[2]
            header = header1 + "|" + header2
            seq_dict[header] = ''
        else:
            seq_dict[header] = seq_dict[header] + line

with open("fos_AA_nogaps_clean.fasta", "w") as fosAA_output_file:
    for header, seq in seq_dict.items():
        fosAA_output_file.write(">"+header+"\n"+seq+"\n")