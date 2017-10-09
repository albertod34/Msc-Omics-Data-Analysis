number_seqs = int(input("How many sequences?: "))
n = 0
seqs_list = []
ind_seqs_list = []

while n < number_seqs:
	seqs = input("Enter the sequence: ")
	ind = seqs.find('ATG')
	seqs_list.append(seqs)
	ind_seqs_list.append(ind)
	n+=1

print(seqs_list)
print(ind_seqs_list)
