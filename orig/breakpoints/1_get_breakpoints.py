####

#is_gtfuse_done  distance        priority        num_cluster     index   this_sample_id  sample  cluster_id      num_disc_read   ch1     start1  stop1   ch2     start2  stop2   orie
ntation     read_orientation        Numts   distance        is_numt_involved        normal_similar  normal_similar_gt3      is_normal_overlapping   tumor_similar   tumor_similar_gt3
      is_tumor_overlap
#1       168     2       5       1       9c274536-3ca1-4f0e-93a8-1688074d862f    BLCA:9c274536-3ca1-4f0e-93a8-1688074d862f_tumor cluster19.1     9       MT      16219   16463   3   
    161285795       161285963       0,1     [9]     HSA_NumtS_163_b1_right  8368694 numt_unlikely   0       0       NA      0       0       t_NA
#1       276     2       5       2       9c274536-3ca1-4f0e-93a8-1688074d862f    BLCA:9c274536-3ca1-4f0e-93a8-1688074d862f_tumor cluster47.1     9       MT      16113   16312   10  
    131045  131321  1,0     [9]     HSA_NumtS_370_b1_left   46711752        numt_unlikely   0       0       NA      0       0       t_NA
#1       275     1       5       3       9c274536-3ca1-4f0e-93a8-1688074d862f    BLCA:9c274536-3ca1-4f0e-93a8-1688074d862f_tumor cluster55.1     49      MT      11174   11448   11  
    67592479        67592754        0,0/0,1 [36, 13]        HSA_NumtS_400_b1_right  9714548 numt_unlikely   0       0       NA      0       0       t_NA
#1       270     1       5       4       9c274536-3ca1-4f0e-93a8-1688074d862f    BLCA:9c274536-3ca1-4f0e-93a8-1688074d862f_tumor cluster56.1     53      MT      11139   11438   11  
    67593972        67594242        1,1     [53]    NA      1000000000      numt_unlikely   0       0       NA      0       0       t_NA
#1       391     2       5       5       9c274536-3ca1-4f0e-93a8-1688074d862f    BLCA:9c274536-3ca1-4f0e-93a8-1688074d862f_tumor cluster63.2     6       MT      16115   16474   15  
    93433252        93433643        1,1/0,0 [1, 5]  NA      1000000000      numt_unlikely   0       0       NA      0       0       t_NA
import sys, subprocess

def CIGARlength(a1):
	a1=a1.replace("N","D").replace("P","D")
	this_countD=a1.count("D")

	while this_countD > 0:
		D_pos=a1.index("D")
		try:
			this_lenD=int(a1[max(0,(D_pos-2)):D_pos])
			a1=a1[0:max(0,(D_pos-2))]+a1[D_pos:len(a1)]
		except:
			a1=a1[0:max(0,(D_pos-1))]+a1[D_pos:len(a1)]
		this_countD=a1.count("D")
	### no deletion now!
	a1=a1.replace("H","M").replace("S","M").replace("I","M")
	a1=a1.split("M")
	this_length=0
	for this_a1 in a1:
		this_length+=int(this_a1)
	return this_length

def CIGAR_removeD(a1):
	a1=a1.replace("N","D").replace("P","D")
	this_countD=a1.count("D")
	while this_countD > 0:
		D_pos=a1.index("D")
		try:
			this_lenD=int(a1[max(0,(D_pos-2)):D_pos])
			a1=a1[0:max(0,(D_pos-2))]+a1[D_pos:len(a1)]
		except:
			a1=a1[0:max(0,(D_pos-1))]+a1[D_pos:len(a1)]
	a1.replace("H","S")
	return a1

def CIGAR_order(a1):
	S_pos=a1.index("S")
	M_pos=a1.index("M")
	if S_pos > M_pos:  #78M23S
		return 1
	else:
		return 0
def CIGAR_Slength(a1):
	Sindex=a1.index("S")
	try:
		Slen=int(a1[(Sindex-2):Sindex)])
	except:
		Slen=int(a1[(Sindex-1):Sindex)])
	return Slen

		


fn="0000_p1p2_candidates.txt"
inputfile=file(fn)
line=inputfile.readline() # header
line=inpurfile.readline()

ofn="r01_"+fn.replace(".txt","_brk.txt")
ofn2="r01_problematic_CIGAR.err"

outputfile=file(ofn,"w")
outputfile2=file(ofn2,"w")

discordant_list=[]

while line:
	outputfile2.write(">>>"+line)

	line_split=line.rstrip().split("\t")
	num_support_reads=int(line_split[8])
	tissue=line_split[6].split(":")[0]
	sample_id=line_split[5]
	cluster_id=line_split[7]
	this_index=int(line_split[4])
	pair_ch=line_split[12]
	pair_start=int(line_split[13])
	pair_stop=int(line_split[14])

	mt_start=int(line_split[10])
	mt_stop=int(line_split[11])

	origin="02_SantaCruz"
	if index >= 100000:
		origin="03_train3"

	### mtDNA translocation file
	mtDNAfn="/nfs/users/nfs_y/ysj/lustreScratch112/MT_pancan/translocation/"+origin+"/"+tissue+"/"+sample_id+"/r02_"+sample_id+"_tumor_MT.MTtr.cluster.txt"
	mtDNAfile=file(mtDNAfn)
	mtDNAline=mtDNAfile.readline()

	num_disc_reads_00=0   ###FLAG
	num_split_reads_00=0
	num_disc_reads_01=0
	num_split_reads_01=0
	num_disc_reads_10=0
	num_split_reads_10=0
	num_disc_reads_11=0
	num_split_reads_11=0

	while mtDNAline:
		mtDNAline_split=mtDNAline.rstrip().split("\t")
		this_cluster_id=mtDNAline_split[-1]
		this_directions=mtDNAline_split[-2]
		this_direction=mtDNAline_split[-2].split(",")[0] # 0 positive, 1 RC
		pair_direction=mtDNAline_split[-2].split(",")[1] # same
		this_CIGAR=mtDNAline_split[5]

		if this_cluster_id != cluster_id:
			mtDNAline=mtDNAfile.readline()
			continue
		discordant_list.append(mtDNAline_split[0:-2])
		
		if ("S" not in this_CIGAR) and ("H" not in this_CIGAR):
			if mtDNAline_split[-2]=="0,0":
				num_disc_reads_00+=1
			elif mtDNAline_split[-2]=="0,1":
				num_disc_reads_01+=1
			elif mtDNAline_split[-2]=="1,0":
				num_disc_reads_10+=1
			elif mtDNAline_split[-2]=="1,1":
				num_disc_reads_11+=1
			else:
				print "Error!Direction?"
				print	mtDNAline
				print this_directions
				raw_input()

		num_clip = this_CIGAR.count("S")+this_CIGAR.count("H")


		if num_clip >=2 : # not suitable for simple breakpoint calculation 
			is_for_brk=0
		else:
			is_for_brk=1

		try:
			this_saz_list=mtDNAline.rstrip().split("\tSA:Z:")[1].split("\t")[0].split(";")[0:-1]
		except:
			mtDNAline=mtDNAfile.readline()
			continue
		this_saz=[]
###### CIGAR
		seq_length=CIGARlength(this_CIGAR)  ### READ LENGTH
		this_CIGAR=CIGAR_removeD(this_CIGAR)  #101M
		this_S_order=CIGAR_order(this_CIGAR)
		this_Slength=CIGAR_Slength(this_CIGAR)
##### FIND breakpoints

		for saz_info in this_saz_list:
			saz_info_split=saz_info.split(",")
			saz_info_ch=saz_info_split[0]
			saz_info_pos=int(saz_info_split[1])
			if saz_info_ch!=pair_ch:
				continue
			if saz_info_pos < pair_start-3000 or saz_info_pos > pair_stop + 3000 ; # proper breakpoint unlikely
				continue
			this_saz.append(saz_info)

		if len(this_saz) == 0:
			mtDNAline=mtDNAfile.readline()
			continue
		if mtDNAline_split[-2]=="0,0":
			num_split_reads_00+=1
		elif mtDNAline_split[-2]=="0,1":
			num_split_reads_01+=1
		elif mtDNAline_split[-2]=="1,0":
			num_split_reads_10+=1
		elif mtDNAline_split[-2]=="1,1":
			num_split_reads_11+=1
		else:
			print "Error!Direction?"
			print	mtDNAline
			print this_directions
			raw_input()

		mtDNA_aligned_pos=int(mtDNAline_split[3])

		try:
			for CIGAR_forM in this_CIGAR_forM:
				if "S" in CIGAR_forM:
					CIGAR_forM=CIGAR_forM.split("S")[1]
				#count
				while 1:
					try:
						is_D_present=CIGAR_forM.index("D")
					except:
						is_D_present=-1
					try:
						is_I_present=CIGAR_forM.index("I")
					except:
						is_I_present=-1
					if is_D_present==-1 and is_I_present==-1:
						break
					if is_D_present > 0 and is_I_present >0: #both present
						if is_D_present < is_I_present:
							len_deleted_bases_mt+=int(CIGAR_forM[0:is_D_present])	
							CIGAR_forM=CIGAR_forM.split("D")[1]
							continue
						else:
							CIGAR_forM=CIGAR_forM.split("I")[1]
							continue
					elif is_D_present > 0:
						len_deleted_bases_mt+=int(CIGAR_forM[0:is_D_present])	
						CIGAR_forM=CIGAR_forM.split("D")[1]
						continue
					else:
						CIGAR_forM=CIGAR_forM.split("I")[1]
						continue
				matched_number_mt+=len(CIGAR_forM)	
		except:
			outputfile.write(mtDNAline)
			mtDNAline=mtDNAfile.readline()
			continue

		mtDNA_brk=mtDNA_aligned_pos+matched_number_mt+len_deleted_bases -1 



	
#	temp_sam="temp_sam.sam"
#	### nuclearDNA translocation file
#	nuDNAfn="/nfs/users/nfs_y/ysj/lustreScratch112/MT_pancan/translocation/05_gtfuse_validation/02_sliced_bam_all/"+str(this_index)+"_nuc_T.bam"
#	subprocess.call("samtools view "+nuDNAfn +" > temp_sam.sam",shell=True)

#	temp_sam_file=file(temp_sam)
#	temp_line=temp_sam_file.readline()

#	while temp_line:
#		if temp_line[0]=="@" : #header
#			temp_line=temp_sam_file.readline()
#			continue	

