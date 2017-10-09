seq = input("Type a sequence: ")

if seq.find("D") == True OR seq.find("E") == True OR seq.find("F") == True OR seq.find("H") == True OR seq.find("I") == True OR seq.find("K") == True OR seq.find("L") == True OR seq.find("M") == True OR seq.find("N") == True OR seq.find("P") == True OR seq.find("Q") == True OR seq.find("R") == True OR seq.find("S") == True OR seq.find("V") == True OR seq.find("W") == True OR seq.find("Y") == True:
    print(seq,"Protein")
elif seq.count("U") >= 1:
    print(seq,"RNA")
else: 
    print(seq, "DNA")
