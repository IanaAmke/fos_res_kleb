from Bio import SeqIO
import time
from collections import defaultdict

start = time.time() 

seq_counts = defaultdict(0)
records = []

for record in SeqIO.parse("fos_AA_final.fasta", "fasta"):  
    seq = str(record.seq)
    if seq in seq_counts:
        records.append(record)
    
    seq_counts[seq] += 1

#writing to a fasta file
SeqIO.write(records, "fos_AA_final_nodup.fasta", "fasta")
end = time.time()

# print the counts
for seq in seq_counts:
    print(seq_counts[seq], seq)

print(f"Run time is {(end- start)/60}") 