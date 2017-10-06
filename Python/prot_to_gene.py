from Bio import ExPASy
from Bio import SwissProt

def prot2gene(prot_name):
	prot_E = ExPASy.get_sprot_raw(prot_name)
	prot_S = SwissProt.read(prot_E)
	print(prot_S.gene_name)
	
	
with open('entries.txt','r') as entries:
	entry = entries.read().split()
counter = 0

for i in entry:
    prot2gene(i)
    counter+=1

	
	

	

	
