import sys
fn1=sys.argv[1] # JuRASS data
fn2=sys.argv[2] # BRASS2 data

ofn="r11_"+"_".join(fn1.split("_")[1:3])+".comparison.txt"
outputfile=file(ofn,"w")


inputfile1=file(fn1)
inputfile2=file(fn2)

lines1=inputfile1.readlines()
lines2=inputfile2.readlines()

len_lines1=len(lines1)
len_lines2=len(lines2)


lines1_split=[]
lines2_split=[]
for i in range(0,len_lines1):
	lines1_this=lines1[i].rstrip().split("\t")
	if lines1_this[3]=="Y" or lines1_this[3]=="MT":
		continue
	if lines1_this[6]=="Y" or lines1_this[6]=="MT":
		continue
	if lines1_this[3]==lines1_this[6]:
		continue
#	lines1_this[4]=int(lines1_this[4])
#	lines1_this[5]=int(lines1_this[5])
#	lines1_this[7]=int(lines1_this[7])
#	lines1_this[8]=int(lines1_this[8])
	lines1_split.append(lines1_this)

for i in range(0,len_lines2):
	if "##" in lines2[i]:
		continue
	lines2_this=lines2[i].rstrip().split("\t")
	if lines2_this[2]=="Y" or lines2_this[2]=="MT":
		continue
	if lines2_this[6]=="Y" or lines2_this[6]=="MT":
		continue
	if lines2_this[2]==lines2_this[6]:
		continue
	lines2_this[4]=int(lines2_this[4])
	lines2_this[5]=int(lines2_this[5])
	lines2_this[8]=int(lines2_this[8])
	lines2_this[9]=int(lines2_this[9])
	lines2_split.append(lines2_this)

len_lines1=len(lines1_split) #total # elements of JuRASS2
len_lines2=len(lines2_split) #total # elements of BRASS2

#print len_lines1
#print len_lines2

J2_only=[]
B2_only=[]
JB_common=[]

for j in range(0,len_lines1):
	line1_this=lines1_split[j]
	line1_brks=line1_this[-2].replace("breakpoints:","").split(";")
	is_J_match=0
	for line1_brk in line1_brks:
		if line1_brk=="":
			continue
		line1_brk_ch1=line1_brk.split(",")[0]
		line1_brk_str1=line1_brk.split(",")[1]
		line1_brk_pos1=int(line1_brk.split(",")[2].split("-(")[0])
		line1_brk_ch2=line1_brk.split(",")[2].split(")-")[1]
		line1_brk_str2=line1_brk.split(",")[3]
		line1_brk_pos2=int(line1_brk.split(",")[4].split("(")[0])

		line1_brk_ch1_revcomp=line1_brk.split(",")[2].split(")-")[1]
		line1_brk_str1_revcomp=line1_brk.split(",")[3].replace("+","_").replace("-","+").replace("_","-")
		line1_brk_pos1_revcomp=int(line1_brk.split(",")[4].split("(")[0])
		line1_brk_ch2_revcomp=line1_brk.split(",")[0]
		line1_brk_str2_revcomp=line1_brk.split(",")[1].replace("+","_").replace("-","+").replace("_","-")
		line1_brk_pos2_revcomp=int(line1_brk.split(",")[2].split("-(")[0])

		for b in range(0,len(lines2_split)):
			line2_this=lines2_split[b]
			line2_brk_ch1=line2_this[2]
			line2_brk_str1=line2_this[3]
			line2_brk_start1=int(line2_this[4])
			line2_brk_stop1=int(line2_this[5])

			line2_brk_ch2=line2_this[6]
			line2_brk_str2=line2_this[7]
			line2_brk_start2=int(line2_this[8])
			line2_brk_stop2=int(line2_this[9])

			if line1_brk_ch1==line2_brk_ch1 and line1_brk_ch2==line2_brk_ch2:
				if not(line1_brk_str1==line2_brk_str1 and line1_brk_str2==line2_brk_str2):
					continue
				if  ((line1_brk_pos1 -  (min(line2_brk_start1,line2_brk_stop1) - 30) ) * (line1_brk_pos1-( max(line2_brk_start1, line2_brk_stop1)+30)) < 0 ) and  ((line1_brk_pos2 -  (min(line2_brk_start2,line2_brk_stop2) - 30) ) * (line1_brk_pos2-( max(lin
e2_brk_start2, line2_brk_stop2)+30)) < 0 )  :   #30bp range, overlap
					is_J_match=1
					break
			elif line1_brk_ch1_revcomp==line2_brk_ch1 and line1_brk_ch2_revcomp==line2_brk_ch2:
				if not(line1_brk_str1_revcomp==line2_brk_str1 and line1_brk_str2_revcomp==line2_brk_str2):
					continue
				if  ((line1_brk_pos1_revcomp -  (min(line2_brk_start1,line2_brk_stop1) - 30) ) * (line1_brk_pos1_revcomp-( max(line2_brk_start1, line2_brk_stop1)+30)) < 0 ) and  ((line1_brk_pos2_revcomp -  (min(line2_brk_start2,line2_brk_stop2) - 30) ) * (
line1_brk_pos2_revcomp-( max(line2_brk_start2, line2_brk_stop2)+30)) < 0 )  :   #30bp range, overlap
					is_J_match=1
					break
			else: # not match
				continue

		if is_J_match==1:  # MATCHED
			JB_common.append(line1_this+["JB_Common:"+line2_this[-1]])
			del lines2_split[b]
			break
	if is_J_match==0: # Not matched
		J2_only.append(line1_this+["Ju_Only"])

		

len_J2_only=len(J2_only)
len_common=len(JB_common)
len_B2_only=len(lines2_split)

for b in range(0,len_B2_only):
	lines2_split[b][4]=str(lines2_split[b][4])
	lines2_split[b][5]=str(lines2_split[b][5])
	lines2_split[b][8]=str(lines2_split[b][8])
	lines2_split[b][9]=str(lines2_split[b][9])
	B2_only.append(lines2_split[b]+["BRASS2_Only"])

outputfile.write("###Results:\tCommon:\t"+str(len_common)+"\t"+";JuOnly:\t"+str(len_J2_only)+"\t"+";BRASS2Only:\t"+str(len_B2_only)+"\n")

for i in range(0,len_common):
	outputfile.write("\t".join(JB_common[i])+"\n")

for i in range(0,len_J2_only):
	outputfile.write("\t".join(J2_only[i])+"\n")

for i in range(0,len_B2_only):
	outputfile.write("\t".join(B2_only[i])+"\n")

print ofn
print("###Results:\tCommon:\t"+str(len_common)+"\t"+";JuOnly:\t"+str(len_J2_only)+"\t"+";BRASS2Only:\t"+str(len_B2_only))
print "done"





