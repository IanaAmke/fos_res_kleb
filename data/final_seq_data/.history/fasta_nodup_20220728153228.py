from Bio import SeqIO
import time

start = time.time() 

seen = set()
records = []

for record in SeqIO.parse("fos_AA_final.fasta", "fasta"):  
    if record.seq not in seen:
        seen.add(record.seq)
        records.append(record)

#writing to a fasta file
SeqIO.write(records, "fos_AA_final_nodup.fasta", "fasta")
end = time.time()

print(f"Run time is {(end- start)/60}") 