with open('input_seqs.txt') as dnas:
	data = dnas.read().splitlines()

n = 0
seqs_list = []
ind_seqs_list = []

while n < len(data):
	ind = data[n].find('ATG')
	seqs_list.append(data)
	ind_seqs_list.append(ind)
	n+=1

print(ind_seqs_list)

ans = ''
while ans != 'Yes' and ans != 'No':
	ans = input("Do you want to see the sequences list?: ")
	if ans == 'Yes':
		print(seqs_list)
	elif ans == 'No':
		print("Bye")
	else:
		print("Please, enter Yes or No")
