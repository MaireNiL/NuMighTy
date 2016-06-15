import sys
fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline()
this_sample=fn.split("_")[1]
ofn="r01_"+"_".join(fn.split("_")[0:2])+".BRASS2num_gt1MB.txt"
outputfile=file(ofn,"w")

count_TR=0
count_IN=0
count_TDP=0
count_TDN=0
count_LD=0

while line:
	if "#EOF" in line:
		break
	line_split=line.rstrip().split("\t")
	ch1=line_split[0]
	ch2=line_split[4]
	pos1=int(line_split[3])
	pos2=int(line_split[7])
	str1=line_split[1]
	str2=line_split[5]
	sample=line_split[12]

	if "," in sample: #germline
		line=inputfile.readline()
		continue

	this_type="NA"

	if ch1!=ch2:
		this_type="translocation"
		count_TR+=1
	else:
		if abs(pos1-pos2) < 1000000:
			line=inputfile.readline()
			continue

		if str1!=str2:
			this_type="inversion"
			count_IN+=1
		elif str1=="+":
			if pos1 > pos2:
				this_type="tandemDuplication_PP"
				count_TDP+=1
			else:
				count_LD+=1
				this_type="largeDeletion_"+str(pos2-pos1+1)
		elif str1=="-":
			if pos2 > pos1:
				this_type="tandemDuplication_NN"
				count_TDN+=1
			else:
				this_type="largeDeletion_"+str(pos1-pos2+1)
				count_LD+=1
	outputfile.write(sample+"\t"+this_type+"\t"+line)
	line=inputfile.readline()

outputfile.write("##\tsample\ttotal\ttranslocation\tinversion\ttandemDup\ttandemDup_PP\ttandemDup_NN\tLargeDeletion\n")
outputfile.write("##\t"+this_sample+"\t"+str(count_TR+count_IN+count_TDP+count_TDN+count_LD)+"\t"+str(count_TR)+"\t"+str(count_IN)+"\t"+str(count_TDP+count_TDN)+"\t"+str(count_TDP)+
"\t"+str(count_TDN)+"\t"+str(count_LD)+"\n")
#outputfile.write("#1_Total\t"+str(count_TR+count_IN+count_TDP+count_TDN+count_LD)+"\t"+sample+"\n")
#outputfile.write("#2_Translocation\t"+str(count_TR)+"\t"+sample+"\n")
#outputfile.write("#2_Inversion\t"+str(count_IN)+"\t"+sample+"\n")
#outputfile.write("#2_TandemDuplidation\t"+str(count_TDP+count_TDN)+"\t"+sample+"\n")
#outputfile.write("#3_TandemDuplidation_PP\t"+str(count_TDP)+"\t"+sample+"\n")
#outputfile.write("#3_TandemDuplidation_NN\t"+str(count_TDN)+"\t"+sample+"\n")
#outputfile.write("#2_LargeDeletion\t"+str(count_LD)+"\t"+sample+"\n")

outputfile.close()
