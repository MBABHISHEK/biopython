## IT is feteching the data from Enterz the gobal biological data base
from Bio import Entrez, SeqIO
Entrez.email = "" 
## Entrez fetch is full form efetch is one of the functions used for reading data and returning as object 
handle = Entrez.efetch(db="nucleotide", id="MN908947", rettype="gb", retmode="text")
##  gene bank format is converted to list format
recs = list(SeqIO.parse(handle, 'gb'))

handle.close()
print(recs)

covid_dna=recs[0].seq
print(covid_dna)

print(f'The genome of Covid-19 consists of {len(covid_dna)} nucleotides.')

# molecular weight
from Bio.SeqUtils import molecular_weight
print("THE  MOLECULAR WEIGHT OF THE COVID DNA IS ",molecular_weight(covid_dna))

#calulate the GUANINE AND CYTOSINE OF TH GIVEN COIVD DNA 
from Bio.SeqUtils import gc_fraction
print("THE  GC CONTENT OF THE COVID DNA IS ", gc_fraction(covid_dna))

#to calculate the number of nucleotides present the given COVID DNA
count_nucleotides = {
    'A': covid_dna.count('A'),
    'T': covid_dna.count('T'),
    'C': covid_dna.count('C'),
    'G': covid_dna.count('G')
}
print(count_nucleotides)

## ploting the graph for the nucleotides

import matplotlib.pyplot as plt
print('graph')
width = 0.5
plt.bar(count_nucleotides.keys(), count_nucleotides.values(), width, color=['b', 'r', 'm', 'c'])
plt.xlabel('Nucleotide')
plt.ylabel('Frequency')
plt.title('Nucleotide Frequency')
plt.show()

covid_mrna = covid_dna.transcribe()
print(covid_mrna)


covid_aa = covid_mrna.translate()
print(covid_aa)



from collections import Counter
common_amino = Counter(covid_aa)

print(common_amino.most_common(10))



del common_amino['*']

width = 0.5
plt.bar(common_amino.keys(), common_amino.values(), width, color=['b', 'r', 'm', 'c'])
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')
plt.title('Protein Sequence Frequency')
plt.show()



print(f"Covid-19's genome has {sum(common_amino.values())} amino acids")


proteins = covid_aa.split('*')

#it's worth to mention that not all the amino acids sequences are proteins. 
#Only the sequences with more than 20 amino acids code for functional proteins.
#The short amino acid sequences are oligopeptides and have other functionalities.
#Here, we will focus on the chains with more than 20 amino acid chains: Proteins.
for protein in proteins:
    if len(protein) < 20:
        proteins.remove(protein)

print(f'We have {len(proteins)} proteins with more than 20 amino acids in the covid-19 genome')
top_5_proteins = sorted(proteins, key = len)

print(top_5_proteins[-1])

print(len(top_5_proteins[-1]))

with open("protein_seq.fasta", "w") as file:
    file.write(f">covid protein\n{top_5_proteins[-1]}")