from Bio.Seq import Seq
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
print("The RNA is:"+"\n"+messenger_rna)
protien=messenger_rna.translate()
print("THE PROTIEN IS  : ",protien)

#the dna to protien
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print("The DNA is:"+"\n"+coding_dna)
protien=coding_dna.translate()
print("THE PROTIEN IS :",protien)