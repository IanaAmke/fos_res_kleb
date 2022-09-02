from email.quoprimime import header_decode
import re

seq_dict = {}

with open("/home/iamke/fos.res.kleb/fos_repo/fos_res_kleb/data/final_seq_data/fos_genes_final.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'fos[A-Z,0-9]*', line)[0]
            print(header1)