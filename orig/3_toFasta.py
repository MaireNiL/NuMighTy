import sys
fn=sys.argv[1]
	
inputfile=file(fn)
line=inputfile.readline()

ofn=fn+".fa"
outputfile=file(ofn,"w")
while line:
	line_split=line.split("\t")
	read_name=line_split[0]+"_flag_"+line_split[1]
	sequence=line_split[9]

	outputfile.write(">"+read_name+"\n"+sequence+"\n")

	line=inputfile.readline()
