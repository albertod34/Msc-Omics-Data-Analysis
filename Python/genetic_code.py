import sequence

seq = input("Enter a DNA, RNA or PROTEIN sequence: ")

print("The sequence is: ", sequence.which_seqtype(seq))


frames=''

if sequence.which_seqtype(seq)== "DNA":
	frames0 = sequence.reading_frames(sequence.dna2rna(seq),0,frames)
	print(sequence.rna2prot(frames0,"gencode.txt"))
	frames1 = sequence.reading_frames(sequence.dna2rna(seq),1,frames)
	print(sequence.rna2prot(frames1,"gencode.txt"))
	frames2 = sequence.reading_frames(sequence.dna2rna(seq),2, frames)
	print(sequence.rna2prot(frames2,'gencode.txt'))
	
	
	
elif sequence.which_seqtype(seq)=="RNA":
	frames0 = sequence.reading_frames(seq,0,frames)
	print(sequence.rna2prot(frames0,"gencode.txt"))
	frames1 = sequence.reading_frames(seq,1,frames)
	print(sequence.rna2prot(frames1,"gencode.txt"))
	frames2 = sequence.reading_frames(seq,2,frames)
	print(sequence.rna2prot(frames2,"gencode.txt"))
	
else:
	print("Proteins have no transformation")
 
