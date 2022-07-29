from difflib import Differ
  
with open('trial.fasta') as file_1, open('trial_2.fasta') as file_2:
    differ = Differ()
  
    for line in differ.compare(file_1.readlines(), file_2.readlines()):
        print(line)