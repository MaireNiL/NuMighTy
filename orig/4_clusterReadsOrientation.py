import sys

def geneinfo(ch1,pos1):
	chr_id=int(ch1.replace("X","23").replace("Y","24").replace("MT","25"))-1
#	print chr_id
	cursor=gene_num[chr_id]/2 # for loose binary searching # 1Mb
	upbound=0
	lowbound=gene_num[chr_id]
	while 1:
		gene_start=infolines[chr_id][cursor][4]
		D=pos1-gene_start
#		print (upbound,lowbound)
		if D >= 0: # read in downstream
			upbound=cursor+1
			cursor=(upbound+lowbound)/2
		elif D < 0 : # read in upstream
			lowbound=cursor
			cursor=(upbound+lowbound)/2
		if upbound >= lowbound:
			break

	this_gene_info=[]

	for ii in range(lowbound,-1,-1):
		try:
			gene_start=infolines[chr_id][ii][4]
		except:
			continue
		gene_stop =infolines[chr_id][ii][5]
		if pos1 < gene_start:
			continue
		if pos1> gene_stop+1000000:
			break
		
		if pos1 < gene_stop:
			this_gene_id=infolines[chr_id][ii][1]
			this_gene_cds_start=infolines[chr_id][ii][6]
			this_gene_cds_stop =infolines[chr_id][ii][7]
			this_gene_exon_num =int(infolines[chr_id][ii][8])
			this_gene_exon_starts =infolines[chr_id][ii][9].split(",")[0:-1]
			this_gene_exon_stops  =infolines[chr_id][ii][10].split(",")[0:-1]
			this_gene_strand=infolines[chr_id][ii][3]
			this_gene_name=infolines[chr_id][ii][12]

			for exonnum in range(0,this_gene_exon_num):
				if pos1 < int(this_gene_exon_starts[exonnum]):
					break
			if this_gene_strand=="-":
				exonnum=this_gene_exon_num-exonnum
			this_gene_info.append([this_gene_name,this_gene_id,str(exonnum),this_gene_strand])
	this_gene_info1_name=[]
	this_gene_info1_all=[]
	for neargene in this_gene_info:
		if neargene[0] not in this_gene_info1_name:
			this_gene_info1_all.append(neargene[0]+"("+neargene[1]+","+neargene[3]+",exonintron"+neargene[2]+")")
			this_gene_info1_name.append(neargene[0])
	return this_gene_info1_all

def output_clusters(cluster,cluster_number):
	cluster1=cluster[:]
	is_generate_output=0
	len_cluster=len(cluster1)
	if len_cluster <= 1:
		return 0
	cluster1.sort(mycmp)
	this_status=0

	subcluster_number=1
	for ii in range(0,len_cluster-1):
		if cluster1[ii][2]==cluster1[ii+1][2] and cluster1[ii][-1]==cluster1[ii+1][-1]:
			if int(cluster1[ii+1][3])-int(cluster1[ii][3]) < 400:
				outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
				this_status=1
				is_generate_output=1
			else:
				if this_status==1:
					outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
					this_status=0
					subcluster_number+=1
		else:
			if this_status==1:
				outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
				this_status=0
				subcluster_number+=1
	if this_status==1: #last line		
		outputfile.write("\t".join(cluster1[len_cluster-1])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
		this_status=0
	return is_generate_output

def mycmp(a1,a2):
	mycmp_a1_ch2=a1[2]
	mycmp_a2_ch2=a2[2]
	if mycmp_a1_ch2=="=":
		mycmp_a1_ch2=a1[2]
	if mycmp_a2_ch2=="=":
		mycmp_a2_ch2=a2[2]
	mycmp_a1_ch2=int(mycmp_a1_ch2.replace("X","23").replace("Y","24").replace("MT","25"))
	mycmp_a2_ch2=int(mycmp_a2_ch2.replace("X","23").replace("Y","24").replace("MT","25"))
	mycmp_a1_orientation=a1[-1]
	mycmp_a2_orientation=a2[-1]
	if mycmp_a1_orientation!=mycmp_a2_orientation:
		return cmp(mycmp_a1_orientation,mycmp_a2_orientation)

	if mycmp_a1_ch2!=mycmp_a2_ch2:
		return cmp(mycmp_a1_ch2,mycmp_a2_ch2)
	mycmp_a1_pos2=int(a1[3])
	mycmp_a2_pos2=int(a2[3])
	return cmp(mycmp_a1_pos2,mycmp_a2_pos2)

#gene_num=[\
#4301,2641,2328,1634,1795,2121,1963,1470,1622,1817,\
#2558,2150,723 ,1371,1419,1655,2391, 622,2774,1180,\
# 544, 924,2117, 348]

#infofnlist=[]
#for i in range(1,26):
#	infofnlist.append("/nfs/users/nfs_y/ysj/Databases/02_genes/01_human/refGene_chr"+str(i).replace("23","X").replace("24","Y").replace("25","MT")+"s.txt")
#infolines=[]
#gene_num=[]
#for i in range(0,25):
#	infolines.append(file(infofnlist[i]).readlines())
#	gene_num.append(len(infolines[i])) #line numbers for gene information
#print gene_num

#for i in range(0,25):
#	for j in range(0,gene_num[i]):
#		infolines[i][j]=infolines[i][j].rstrip().split("\t")
#		infolines[i][j][4]=int(infolines[i][j][4])
#		infolines[i][j][5]=int(infolines[i][j][5])
fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline()

ofn="./00006_read_clusters/r06_"+fn.replace(".txt",".clusters_ori.txt")
outputfile=file(ofn,"w")


while line[0] == "@": #header
	if line[0:3]=="@RG":
		outputfile.write(line) 
	line=inputfile.readline()		
#end of header
#IL16_5041:1:35:8908:5728        65      MT      11360   25      108M    X       5496252 0       ATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTAAACGCT    DBGBGGHHGHDHHHHG>GGGH
HH@HHHHHHHHHHDHHFHHHHHHGHHDHHHHGFFHEHHDECGBGGADD<DDDGGDGBC8=1?88=></=AACAA>BB@B4>>@?/:&    X0:i:1  X1:i:0  MD:Z:102C0T1T0T0A0      RG:Z:1059070    XG:i:0  AM:i:25 NM:i:5  SM:i:25 XM:i:5  XO:i:0  XT:A:U
prev_pos_read1=-1000
#prev_line_status=0 # when clustered 1, otherwise 0
cluster_number=0
this_cluster=[]
while line:
	line_split=line.rstrip().split("\t")
	read_id=line_split[0]
	read_flag=int(line_split[1])
	read2_chr=line_split[2]
	read2_pos=int(line_split[3])

	read1_chr=line_split[6]
	if read1_chr=="=":
		read2_chr=read1_chr
	read1_pos=int(line_split[7])

#	read1_gene=geneinfo(read1_chr,read1_pos)
#	raw_input()
#	if read1_gene==[]:
#		line=inputfile.readline()
#		continue
#	read2_gene=geneinfo(read2_chr,read2_pos)
#	raw_input()
#	if read2_gene==[]:
#		line=inputfile.readline()
#		continue
#	this_line_info=line.replace("\n","\t")+",".join(read1_gene)+"\t"+",".join(read2_gene)+"\n"
#	this_line_info=line_split+[",".join(read1_gene),",".join(read2_gene)]
#	if read1_pos - prev_pos_read1 < 0: # Chr changed
#		if prev_line_status==1:
#			outputfile.write(prev_line_info)
#			prev_line_status=0
#	elif read1_pos - prev_pos_read1 < 400: #clustered
#		outputfile.write(prev_line_info)
#		prev_line_status=1
#	else:
#		if prev_line_status==1:
#			outputfile.write(prev_line_info)
#			prev_line_status=0
#	prev_line_info=this_line_info
#	prev_pos_read1=read1_pos

	is_this_rc=read_flag%256/16%8%4%2
	is_next_rc=read_flag%256/16%8%4/2

#	if "+" in ",".join(read1_gene) and "-" not in ",".join(read1_gene):
#		ori_gene1=0
#	elif "-" in ",".join(read1_gene) and "+" not in ",".join(read1_gene):
#		ori_gene1=1
#	else:
#		ori_gene1=2 #complex
#	if "+" in ",".join(read2_gene) and "-" not in ",".join(read2_gene):
#		ori_gene2=0
#	elif "-" in ",".join(read2_gene) and "+" not in ",".join(read2_gene):
#		ori_gene2=1
#	else:
#		ori_gene2=2 #complex
#	fusion_type="complex_conc" #concordant
#	if max(ori_gene1,ori_gene2)<2:
#		orientations=(ori_gene1,ori_gene2,is_this_rc,is_next_rc)
#		if orientations in [(0,0,0,1),(0,1,0,0),(1,0,1,1),(1,1,1,0)]:
#			fusion_type="f1f2_conc" #intact fusion   ==1==>  ==2==>
#		elif orientations in [(0,0,1,0),(0,1,1,1),(1,0,0,0),(1,1,0,1)]:
#			fusion_type="r1r2_conc" #f2f1 # intact fusion  <==1== <==2==
#		elif orientations in [(0,0,0,0),(0,1,0,1),(1,0,1,0),(1,1,1,1)]:
#			fusion_type="f1r2_disc" #f1r2 # XXX ==1==> <==2==
#		elif orientations in [(0,0,1,1),(0,1,1,0),(1,0,0,1),(1,1,0,0)]:
#			fusion_type="r1f2_disc" #f2f1 # intact fusion  <==1== ==2==>
#		else:
#			print("Orientation error!")
#			raw_input()
	fusion_type=str(is_this_rc)+","+str(is_next_rc)
#	this_line_info=line_split+[",".join(read1_gene),",".join(read2_gene),fusion_type]
	this_line_info=line_split+[fusion_type]
	if read1_pos - prev_pos_read1 > 400 or read1_pos - prev_pos_read1 < 0: # Chr changed or not in the same cluster
		#TREAT THIS CLUSTER
		aa=output_clusters(this_cluster,cluster_number)
		if aa==1:
			cluster_number+=1
		this_cluster=[this_line_info]
	elif read1_pos - prev_pos_read1 <= 400: #clustered
		this_cluster.append(this_line_info)

	prev_pos_read1=read1_pos
#	outputfile.write(line.replace("\n","\t")+",".join(read1_gene)+"\t"+",".join(read2_gene)+"\n")
	line=inputfile.readline()

#TREAT THIS CLUSTER
output_clusters(this_cluster,cluster_number)

