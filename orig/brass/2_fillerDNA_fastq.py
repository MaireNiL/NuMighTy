import sys
fn=sys.argv[1]

inputfile=file(fn)
line=inputfile.readline()

ofn=fn.replace("r01","r02").replace(".txt","_filler.fa")
outputfile=file(ofn,"w")

i=0
while line:
	i+=1
	if "#" in line[0]:
		outputfile.close()
		break
	line_split=line.rstrip().split("\t")
	
	seq=line_split[10]

	seq_id=">"+line_split[0]+"_"+line_split[1]+"_"+str(i)	
	seq_len=len(seq)

	if seq==".":
		seq="NNNNNNNNNN"
		seq_len=0

	outputfile.write(seq_id+"\n"+seq+"\n")
	line=inputfile.readline()


