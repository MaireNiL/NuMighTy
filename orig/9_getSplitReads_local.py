import sys, os

fn1=sys.argv[1]

project=sys.argv[2]

sample_id=fn1.split("_")[2]


inputfile=file(fn1)
line=inputfile.readline()

ofn=fn1.replace("r08_","r09_").replace(".ctr","_brk1.ctr")
outputfile=file(ofn,"w")
distance_thr=500

total_cases=0
brk_cases=0
while line:
	if "match" in line:  # JuRASS 1, not somatic
		line=inputfile.readline()
		continue

	total_cases+=1
	line_split=line.rstrip().split("\t")
	ch1=line_split[3]
	start1=int(line_split[4])
	stop1=int(line_split[5])
	distance1=stop1-start1+1


	ch2=line_split[6]
	start2=int(line_split[7])
	stop2=int(line_split[8])
	distance2=stop2-start2+1

	orientation=line_split[9].split(",")
	print("\n\n*************")
	print line
	if distance1 > 2000 or distance2 > 2000: #not real cluster
		outputfile.write(line.replace("\n","\tnot_a_cluster\n"))
		print("...not a cluster...")
		line=inputfile.readline()
		continue
	
	if orientation[0]=="0":
		l_brk_upper=(stop1+start1)/2
		l_brk_lower=(stop1+distance_thr)
		l_spl_upper=start1-distance_thr
		l_spl_lower=stop1+100
	elif orientation[0]=="1":
		l_brk_lower=(stop1+start1)/2
		l_brk_upper=(start1-distance_thr)
		l_spl_upper=start1-100
		l_spl_lower=stop1+distance_thr
	else:
		print("No orientation!")
		print line
		raw_input()

	if orientation[1]=="0":
		r_brk_upper=(stop2+start2)/2
		r_brk_lower=(stop2+distance_thr)
		r_spl_upper=start2-distance_thr
		r_spl_lower=stop2+100
	elif orientation[1]=="1":
		r_brk_lower=(stop2+start2)/2
		r_brk_upper=(start2-distance_thr)
		r_spl_upper=start2-100
		r_spl_lower=stop2+distance_thr
	else:
		print("No orientation!")
		print line
		raw_input()

	###GET UNMAPPED READS
	print("...get unmmaped reads")
	if orientation[0]=="0":
		command_l="samtools view -f 4 -F 40 "
	else:
		command_l="samtools view -f 36 -F 8 "
	command_l=command_l+" /nfs/cancer_ref01/nst_links/live/"+project+"/"+sample_id+"/"+sample_id+".bam "+ch1+":"+str(l_spl_upper)+"-"+str(l_spl_lower)+ " > "+sample_id+".splits.temp"
	if orientation[1]=="0":
		command_r="samtools view -f 4 -F 40 "
	else:
		command_r="samtools view -f 36 -F 8 "
	command_r=command_r+" /nfs/cancer_ref01/nst_links/live/"+project+"/"+sample_id+"/"+sample_id+".bam "+ch2+":"+str(r_spl_upper)+"-"+str(r_spl_lower)+ " >> "+sample_id+".splits.temp"

	
	print command_l
	os.system(command_l)
	print command_r
	os.system(command_r)

	
	###CONVERT to fasta	
	print("...convert to fasta")
	command_convert="python 03_toFasta.py "+sample_id+".splits.temp"
	os.system(command_convert)

	tempfn=sample_id+".splits.temp"
	tempfile=file(tempfn)
	rc_temp=0
	templine=tempfile.readline()

	while templine:
		rc_temp+=1
		templine=tempfile.readline()

		if rc_temp>3000:
			break
	if rc_temp>3000:
		outputfile.write(line.replace("\n","\tnot_a_cluster(too_many_reads)\n"))
		print("...not a cluster...")
		line=inputfile.readline()
		continue
	if rc_temp==0:
		outputfile.write(line.replace("\n","\tnot_unmapped_reads\n"))
		print("...no unmapped reads...")
		line=inputfile.readline()
		continue

	###PREPARE TEMPLATE
	l_temp_upper=max(start1-2000,1)
	l_temp_lower=stop1+2000  # for BLAT
	r_temp_upper=max(start2-2000,1)
	r_temp_lower=stop2+2000
     
	dbfn=sample_id+".database_blat.temp.fa"
	output_db=file(dbfn,"w")
	output_db.write(">"+ch1+"\n")
 
	tempfn1="/nfs/users/nfs_y/ysj/Databases/01_references/01_hg19/chr"+ch1+".fa"
	tempfile1=file(tempfn1)
	templine1=tempfile1.readline()
 
	offset1=(l_temp_upper-1)/50*51+(l_temp_upper-1)%50
	tempfile1.seek(offset1,1)
 
	for base_number in range(l_temp_upper,l_temp_lower):
	        base=tempfile1.read(1)
	        if base=="\n":
		       output_db.write(base)
		       base=tempfile1.read(1)
	        output_db.write(base)
	output_db.write("\n")
 
	output_db.write(">"+ch2+"\n")
     
	tempfn2="/nfs/users/nfs_y/ysj/Databases/01_references/01_hg19/chr"+ch2+".fa"
	tempfile2=file(tempfn2)
	templine2=tempfile2.readline()
	
	offset2=(r_temp_upper-1)/50*51+(r_temp_upper-1)%50
	tempfile2.seek(offset2,1)
 
	for base_number in range(r_temp_upper,r_temp_lower):
	        base=tempfile2.read(1)
	        if base=="\n":
		       output_db.write(base)
		       base=tempfile2.read(1)
	        output_db.write(base)
	output_db.write("\n")
	output_db.close()


	###BLAT
	print("...blat")
#	dbfn=sample_id+".database_blat.temp"
#	output_db=file(dbfn,"w")
#	output_db.write("/nfs/users/nfs_y/ysj/Databases/01_references/01_hg19/chr"+ch1+".fa\n")
#	output_db.write("/nfs/users/nfs_y/ysj/Databases/01_references/01_hg19/chr"+ch2+".fa\n")

	command_blat="~/Workarea/blat_250613/blat -minScore=10  -tileSize=8 "+sample_id+".database_blat.temp.fa "+sample_id+".splits.temp.fa "+sample_id+".splits.temp.fa.psl"
	os.system(command_blat)

	print("...blat done. Interpreting blat results...")
	###BLAT interpret

	blatfn=sample_id+".splits.temp.fa.psl"
	blatfile=file(blatfn)
	blatline=blatfile.readline() #header
	blatline=blatfile.readline() #header
	blatline=blatfile.readline() #header
	blatline=blatfile.readline() #header
	blatline=blatfile.readline() #header
	blatline=blatfile.readline() #header

	prev_read_name=""
	blatlines1=[]
	blatlines2=[]
	is_ch1_in=0
	is_ch2_in=0
	breakpoints=[]
	breakpoints_frequency=[]
	while blatline:
		blatline_split=blatline.rstrip().split("\t")

		try:
			blatline_name=blatline_split[9]
			blat_str=blatline_split[8]
			query_start=int(blatline_split[11])
			query_stop=int(blatline_split[12])
			template_ch=blatline_split[13].replace("chr","").replace("MT","M").replace("M","MT")
		except:
			blatline=blatfile.readline()
			continue
		if query_start <= 100-query_stop:  #----------| read_length=100
			read_type=0  # left fragment
		elif query_start > 100 -query_stop:  # |----------
			read_type=1 # right fragment
		
		breakpoint=-1
		if blat_str=="+":
			if read_type==0:
				if template_ch==ch1:
					breakpoint=int(blatline_split[16])+l_temp_upper-1
				elif template_ch==ch2:
					breakpoint=int(blatline_split[16])+r_temp_upper-1
				else:
					print("Problem! CH not matched...!")
					raw_input()
			else:
				if template_ch==ch1:
					breakpoint=int(blatline_split[15])+l_temp_upper
				elif template_ch==ch2:
					breakpoint=int(blatline_split[15])+r_temp_upper
				else:
					print("Problem! CH not matched...!")
					raw_input()
		else:
			if read_type==1:
#				breakpoint=int(blatline_split[16])
				if template_ch==ch1:
					breakpoint=int(blatline_split[16])+l_temp_upper-1
				elif template_ch==ch2:
					breakpoint=int(blatline_split[16])+r_temp_upper-1
				else:
					print("Problem! CH not matched...!")
					raw_input()
			else:
				#breakpoint=int(blatline_split[15])+1
				if template_ch==ch1:
					breakpoint=int(blatline_split[15])+l_temp_upper
				elif template_ch==ch2:
					breakpoint=int(blatline_split[15])+r_temp_upper
				else:
					print("Problem! CH not matched...!")
					raw_input()

		if blatline_name==prev_read_name:
			if (ch1==template_ch and breakpoint <= l_brk_lower and breakpoint >= l_brk_upper) :
				is_ch1_in=1
				blatlines1.append([blatline_name,query_start,query_stop,template_ch,breakpoint,read_type,blat_str,int(blatline_split[15]),int(blatline_split[16])])
			elif (ch2==template_ch and breakpoint <= r_brk_lower and breakpoint >= r_brk_upper):
				is_ch2_in=1
				blatlines2.append([blatline_name,query_start,query_stop,template_ch,breakpoint,read_type,blat_str,int(blatline_split[15]),int(blatline_split[16])])
			blatline=blatfile.readline()

		else: ##analyse blatlines
			if is_ch1_in * is_ch2_in == 1:
				if len(blatlines1) == 1 or len(blatlines2) == 1 and blatlines1[0][5]!=blatlines2[0][5]:

					query_len=(blatlines1[0][2]-blatlines1[0][1])+(blatlines2[0][2]-blatlines2[0][1])
					if query_len > 90:
						base_overlap=min(blatlines1[0][2],blatlines2[0][2])-max(blatlines1[0][1],blatlines2[0][1]) # if 0, no microhomology; > 0 n bp homology,; < 0, base insertion
						if base_overlap < 30:
							if blatlines1[0][5]==0: # 1 is left
								this_breakpoint=blatlines1[0][3]+","+blatlines1[0][6]+","+str(blatlines1[0][4])+"-("+str(base_overlap)+")-"+blatlines2[0][3]+","+blatlines2[0][6]+","+str(blatlines2[0][4])
							else: #1 is right
								this_breakpoint=blatlines1[0][3]+","+blatlines1[0][6].replace("+","_").replace("-","+").replace("_","-")+","+str(blatlines1[0][4])+"-("+str(base_overlap)+")-"+blatlines2[0][3]+","+blatlin
es2[0][6].replace("+","_").replace("-","+").replace("_","-")+","+str(blatlines2[0][4])

							if this_breakpoint not in breakpoints:
								breakpoints.append(this_breakpoint)
								breakpoints_frequency.append(1)
							else:
								breakpoints_frequency[breakpoints.index(this_breakpoint)]+=1

			prev_read_name=blatline_name
			blatlines1=[]
			blatlines2=[]
			is_ch1_in=0
			is_ch2_in=0
			if (ch1==template_ch and breakpoint <= l_brk_lower and breakpoint >= l_brk_upper) :
				is_ch1_in=1
				blatlines1.append([blatline_name,query_start,query_stop,template_ch,breakpoint,read_type,blat_str,int(blatline_split[15]),int(blatline_split[16])])
			elif (ch2==template_ch and breakpoint <= r_brk_lower and breakpoint >= r_brk_upper):
				is_ch2_in=1
				blatlines2.append([blatline_name,query_start,query_stop,template_ch,breakpoint,read_type,blat_str,int(blatline_split[15]),int(blatline_split[16])])
			blatline=blatfile.readline()

	#LASTline
	if is_ch1_in * is_ch2_in == 1:
		if len(blatlines1) == 1 or len(blatlines2) == 1 and blatlines1[0][5]!=blatlines2[0][5]:
			query_len=(blatlines1[0][2]-blatlines1[0][1])+(blatlines2[0][2]-blatlines2[0][1])
			if query_len > 90:

				base_overlap=min(blatlines1[0][2],blatlines2[0][2])-max(blatlines1[0][1],blatlines2[0][1]) # if 0, no microhomology; > 0 n bp homology,; < 0, base insertion
				if base_overlap < 30:
					if blatlines1[0][5]==0: # 1 is left
						this_breakpoint=blatlines1[0][3]+","+blatlines1[0][6]+","+str(blatlines1[0][4])+"-("+str(base_overlap)+")-"+blatlines2[0][3]+","+blatlines2[0][6]+","+str(blatlines2[0][4])
					else: #1 is righ
						this_breakpoint=blatlines1[0][3]+","+blatlines1[0][6].replace("+","_").replace("-","+").replace("_","-")+","+str(blatlines1[0][4])+"-("+str(base_overlap)+")-"+blatlines2[0][3]+","+blatlines2[0][6].replac
e("+","_").replace("-","+").replace("_","-")+","+str(blatlines2[0][4])
		#			this_breakpoint=blatlines1[0][3]+","+blatlines1[0][6]+","+str(blatlines1[0][4])+"-("+str(base_overlap)+")-"+blatlines2[0][3]+","+blatlines2[0][6]+","+str(blatlines2[0][4])
					if this_breakpoint not in breakpoints:
						breakpoints.append(this_breakpoint)
						breakpoints_frequency.append(1)
					else:
						breakpoints_frequency[breakpoints.index(this_breakpoint)]+=1


	###OUTPUT RESULTS
	print_contents=""
	if breakpoints_frequency==[]:
		outputfile.write(line.replace("\n","\tnot_a_splitread\n"))
		print(".....***not_a_splitread***"+str(brk_cases)+"/"+str(total_cases))
	else:
		brk_cases+=1
		breakpoints.sort(key=dict(zip(breakpoints,breakpoints_frequency)).get,reverse=True)
		breakpoints_frequency.sort(reverse=True)
		outputfile.write(line.replace("\n","\tbreakpoints:"))
#		print breakpoints
#		print breakpoints_frequency
	
		for n,x in enumerate(breakpoints_frequency):
			outputfile.write(breakpoints[n]+"("+str(x)+");")
			print_contents+=breakpoints[n]+"("+str(x)+");"
		outputfile.write("\n")
		print(".....****BREAKPOINTS***"+print_contents+" "+str(brk_cases)+"/"+str(total_cases))
#	os.system("rm "+sample_id+".*.temp*")

	line=inputfile.readline()

print("Done")
outputfile.close()





ml28@cgpbar:/nfs/users/nfs_y/ysj/lustreScratch112/MT_translocations/06_rearrangements_nucleus/00006_read_clusters$  
