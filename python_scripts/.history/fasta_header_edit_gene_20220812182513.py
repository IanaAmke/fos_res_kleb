import re

seq_dict = {}

with open("/home/iamke/fos.res.kleb/fos_repo/fos_res_kleb/data/final_seq_data/fos_genes_final.fasta") as fosgene_file:
    for line in fosgene_file:
        line = line.rstrip()
        if line.startswith(">"):
            header1 = re.findall(r'fos[A-Z,0-9]*', line)[0]
            header2 = re.split('\s', line)[0]
            if line.startswith(">fos"):
                header2 = re.split('|', line)[1]
            if line.startswith(">gb"):
                header2 = re.split('|', line)[2]
            header = header1 + "|" + header2
            seq_dict[header] = ''
        else:
            seq_dict[header] = seq_dict[header] + line

with open("fos_genes_final_clean.fasta", "w") as fosgene_output_file:
    for header, seq in seq_dict.items():
        fosgene_output_file.write(">"+header+"\n"+seq+"\n")