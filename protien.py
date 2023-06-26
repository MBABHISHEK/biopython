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
width = 0.5
plt.bar(count_nucleotides.keys(), count_nucleotides.values(), width, color=['b', 'r', 'm', 'c'])
plt.xlabel('Nucleotide')
plt.ylabel('Frequency')
plt.title('Nucleotide Frequency')