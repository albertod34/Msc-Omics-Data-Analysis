import sequence
import Bio
from Bio.Data import CodonTable
codons = CodonTable.unambiguous_dna_by_name['Standard']


def prot2dna(seq):
    dna=''
    dic = codons.back_table
    for aas in seq:
        for keys in dic:
            if keys == aas:
                dna += dic[keys]
    return(dna)
                

seq = raw_input("Enter a DNA, RNA or PROTEIN sequence: ") #The argument input gives me problems in python2, so i used raw_input (I couldn't install BioPython in python3)

print("The sequence is: ", sequence.which_seqtype(seq))



if sequence.which_seqtype(seq)== "DNA":
    print(seq)
    print(sequence.dna2rna(seq))
    print(sequence.rna2prot(sequence.dna2rna(seq),'gencode.txt'))
    

    
    
    
elif sequence.which_seqtype(seq)=="RNA":
    print(seq.replace("U","T"))
    print(seq)
    print(sequence.rna2prot(seq,'gencode.txt'))

    
else:
    print(prot2dna(seq))
    print(sequence.dna2rna(prot2dna(seq)))
    print(seq)
    
