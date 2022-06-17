import os,sys,glob

if __name__ == "__main__":

	infile = sys.argv[1]
	temp = []
	outfile = "temp.pdb"
	with open(infile,"r") as F:

		for line in F.readlines():
			tline = line.strip()
			if tline.startswith("ATOM") or tline.startswith("ANISOU"):
				oline = tline[:21] + "A" + tline[22:]
			else:
				oline = tline
			print(oline)
			temp.append(oline)
	print(temp)
	with open(outfile,"w") as W:
		for line in temp:
			W.write(line+"\n")
