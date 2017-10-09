seqs_list = []
seq = 'waiting seqs'

while seq != '':
	seq = input("Enter the sequence: ")
	if 'U' in seq:
		print("Only DNA sequences are allowed")
	else:
		seqs_list.append(seq)

seqs_list.pop()
print(seqs_list)
