#THIS script allows cluster with single supporting read
import sys

def output_clusters(cluster,cluster_number):
	cluster1=cluster[:]
	is_generate_output=0
	len_cluster=len(cluster1)
#	if len_cluster <= 1:    # for normal: output cluster even == 1
	if len_cluster < 1:    # for normal: output cluster even == 1
		return 0
	cluster1.sort(mycmp)
	this_status=0

	subcluster_number=1
	if len_cluster==1:
		outputfile.write("\t".join(cluster1[0])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
		is_generate_output=1
		return is_generate_output


	for ii in range(0,len_cluster-1):
		if cluster1[ii][2]==cluster1[ii+1][2] and cluster1[ii][-1]==cluster1[ii+1][-1]:
			if int(cluster1[ii+1][3])-int(cluster1[ii][3]) < 500:
				outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
				this_status=1
				is_generate_output=1
			else:
#				if this_status==1: #for normal: output cluster even -==1
				if this_status<=1:
					outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
					this_status=0
					subcluster_number+=1
		else:
#			if this_status==1: #for normal: output cluster even -==1
			if this_status<=1:
				outputfile.write("\t".join(cluster1[ii])+"\tcluster"+str(cluster_number+1)+"."+str(subcluster_number)+"\n")
				this_status=0
				subcluster_number+=1
#	if this_status==1: #last line	#for normal: output cluster even ==1
	if this_status<=1:
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

fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline()

ofn=sys.argv[2]
#ofn="./00006_read_clusters/r06_"+fn.replace(".txt",".clusters_ori.txt")
outputfile=file(ofn,"w")


while line[0] == "@": #header
	if line[0:3]=="@RG":
		outputfile.write(line) 
	line=inputfile.readline()		
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

	is_this_rc=read_flag%256/16%8%4%2
	is_next_rc=read_flag%256/16%8%4/2

	fusion_type=str(is_this_rc)+","+str(is_next_rc)
	this_line_info=line_split+[fusion_type]
	if read1_pos - prev_pos_read1 > 500 or read1_pos - prev_pos_read1 < 0: # Chr changed or not in the same cluster
		#TREAT THIS CLUSTER
		aa=output_clusters(this_cluster,cluster_number)
		if aa==1:
			cluster_number+=1
		this_cluster=[this_line_info]
	elif read1_pos - prev_pos_read1 <= 500: #clustered
		this_cluster.append(this_line_info)

	prev_pos_read1=read1_pos
	line=inputfile.readline()

output_clusters(this_cluster,cluster_number)
