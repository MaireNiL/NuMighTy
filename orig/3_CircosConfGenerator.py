
samplefn="3_example_sample_list.txt"
inputfile=file(samplefn)
line=inputfile.readline()

while line:
	line_split=line.rstrip().split("\t")
	cfn=line_split[0]
	nfn=line_split[1]

	this_sample="_".join(cfn.split("_")[0:2])

	template="PD4107.conf"
	tempfile=file(template)
	templine=tempfile.readline()

	ofn=this_sample+".conf"
	outputfile=file(ofn,"w")

	while templine:
		if "###" not in templine:
			outputfile.write(templine)
			templine=tempfile.readline()
			continue
		if "###OUTPUTFILE###" in templine:
			outputfile.write(templine)
			templine=tempfile.readline()
			templine=tempfile.readline()
			outputfile.write("file\t= "+this_sample+".png\n")
			continue
		elif "###TUMOR" in templine:
			outputfile.write(templine)
			templine=tempfile.readline()
			templine=tempfile.readline()
			outputfile.write("file\t= "+cfn+"\n")
			continue
		elif "###NORMAL" in templine:
			outputfile.write(templine)
			templine=tempfile.readline()
			templine=tempfile.readline()
			outputfile.write("file\t= "+nfn+"\n")
			continue

		print "PROBLEM!"
		print templine
		raw_input()

	line=inputfile.readline()
