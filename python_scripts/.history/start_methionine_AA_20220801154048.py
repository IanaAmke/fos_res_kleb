from Bio import SeqIO

filename = "/home/iamke/fos.res.kleb/fos_repo/fos_res_kleb/data/final_seq_data/fos_AA_final.fasta"

bad = 0

for record in SeqIO.parse(filename, "fasta"):
    if not record.seq.startswith("M"):
        bad = bad + 1
        print(record.id + " starts " + record.seq[0])
print("Found " + str(bad) + " records in " + filename + " which did not start with M")