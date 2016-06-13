import sys

fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline() #header

sample_name="_".join(fn.split("_")[1:3])

ofn=fn.replace("r06_","r07_").replace(".txt",".ctr")
outputfile=file(ofn,"w")
while "@RG" in line[0:3]:
	outputfile.write(line) #header
	line=inputfile.readline()
prev_cluster_id=""
element_num=0
read1_chr=[]
read2_chr=[]
read1_start=[]
read2_start=[]
#read1_genes=[]
#read2_genes=[]
orientations=[]
orientation_num=[]
#line=inputfile.readline()

while line:
	line_split=line.rstrip().split("\t")
	cluster_id=line_split[-1]

	if cluster_id!=prev_cluster_id:
		if prev_cluster_id!="":
			outputfile.write(sample_name+"\t"+prev_cluster_id+"\t"+str(element_num)+"\t"+\
read1_chr[0]+"\t"+str(min(read1_start))+"\t"+str(max(read1_start))+"\t"+\
read2_chr[0]+"\t"+str(min(read2_start))+"\t"+str(max(read2_start))+"\t"+"/".join(orientations)+"\t"+str(orientation_num)+"\n")
		prev_cluster_id=cluster_id
		element_num=1
		read1_chr=[line_split[2]]
		read2_chr=[line_split[6].replace("=",line_split[2])]
		read1_start=[int(line_split[3])]
		read2_start=[int(line_split[7])]
#		read1_genes=[line_split[-4]]
#		read2_genes=[line_split[-3]]
		orientations=[line_split[-2]]
		orientation_num=[1]
		line=inputfile.readline()
		continue
	read1_start.append(int(line_split[3]))
	read2_start.append(int(line_split[7]))
#	if line_split[-4] not in read1_genes:
#		read1_genes.append(line_split[-4])
#	if line_split[-3] not in read2_genes:
#		read2_genes.append(line_split[-3])
	if line_split[-2] not in orientations:
		orientations.append(line_split[-2])
		orientation_num.append(1)
	else:
		orientation_num[orientations.index(line_split[-2])]+=1
	element_num+=1
	line=inputfile.readline()
	
outputfile.write(sample_name+"\t"+prev_cluster_id+"\t"+str(element_num)+"\t"+\
read1_chr[0]+"\t"+str(min(read1_start))+"\t"+str(max(read1_start))+"\t"+\
read2_chr[0]+"\t"+str(min(read2_start))+"\t"+str(max(read2_start))+"\t"+"/".join(orientations)+"\t"+str(orientation_num)+"\n")
