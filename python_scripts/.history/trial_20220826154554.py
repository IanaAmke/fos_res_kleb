from email.quoprimime import header_decode
import re

seq_dict = {}

with open("/home/iamke/fos_ncbi/fos_genes_fasta/fos_genes_all_db_NR_nogaps.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'fos[A-Z,0-9]*|Orf[A-Z,0-9]'
            print(header1)