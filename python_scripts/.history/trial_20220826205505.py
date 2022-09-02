from email.quoprimime import header_decode
import re

seq_dict = {}

with open("/home/iamke/fos_ncbi/fos_genes_fasta/fos_genes_all_db_NR_nogaps.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'fos[A-Z,0-9]*|Fos[A-Z,0-9]', line)[0]
            header2 = re.split('\s', line)[1]
            header = header1 + "|" + header2
            seq_dict[header] = ''
        else:
            seq_dict[header] = seq_dict[header] + line

with open("fos_genes_all_db_clean.fasta", "w") as fosgene_output_file:
    for header, seq in seq_dict.items():
        fosgene_output_file.write(">"+header+"\n"+seq+"\n")