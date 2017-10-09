with open('input_seqs.txt') as dnas:
	data = dnas.read().splitlines()

n = 0
seqs_list = []

while n < len(data):
	ind = data[n].find('ATG')
	if ind != -1:
		seqs_list.append(data[n][ind:])
	n+=1

print(seqs_list)
