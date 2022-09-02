import sys
import glob
import csv
import os
from Bio import SeqIO


included_fos=[]
uniquehits_file=[]
for file in glob.glob("fosgeneHits*.fasta"):
	for seq_record in SeqIO.parse(file, "fasta"):
		header = str(seq_record.id)
		sequence = str(seq_record.seq).upper()
		fasta_name=(header.split("/")[6])
		fasta_name=fasta_name.replace("_NextSeq_500_paired_end_sequencing","")
		fasta_name=fasta_name.replace("_Illumina_MiSeq_paired_end_sequencing","")

		fasta1=fasta_name.split(".fasta_")
		with open("Huynh2020_sample_alias.tsv","r") as huynh:
			tsv_file = csv.reader(huynh, delimiter="\t")
			next(tsv_file, None)
			for line in tsv_file:
				if fasta1[0] == line[0]:
					fasta1[0]=line[1]

		with open("fosgenes_antibiogram_uniqueHits.tsv", "r") as antibiogram:
			antibiogram_file = csv.reader(antibiogram, delimiter="\t")
			next(antibiogram_file, None)
			for line in antibiogram_file:
				unique_hit=[]
				unique_hit.append(line[0])
				unique_hit.append(line[15])
				uniquehits_file.append(unique_hit)
				if line[0] == fasta1[0] and line[15] == fasta1[1]:
					if fasta1 not in included_fos:
						with open("filtered_fosgenes.fasta","a") as new_file:
							new_file.write(">" + line[0] + "_" + line[15] + "_"+ line[14] + "\n" + sequence + "\n")	
							included_fos.append(fasta1)		


first_set = set(map(tuple, included_fos))
secnd_set = set(map(tuple, uniquehits_file))

print("These are the missing genes:")
print(first_set.symmetric_difference(secnd_set))

