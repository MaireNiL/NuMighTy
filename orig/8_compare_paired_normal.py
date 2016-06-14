import sys
fn1=sys.argv[1]
fn2=sys.argv[2]

inputfile1=file(fn1)
line1=inputfile1.readline()

inputfile2=file(fn2)
line2=inputfile2.readline()

ofn=fn1.replace("r07_","r08_").replace(".ctr","_pn.ctr")
outputfile=file(ofn,"w")

distance_thr=1000

ii=0
while line1:
	ii+=1
	if ii%100==0:
		print ii
		print line1
	line1_split=line1.rstrip().split("\t")
	ch1_1=int(line1_split[3].replace("X","0").replace("Y","-1").replace("MT","-2"))
	ch1_2=int(line1_split[6].replace("X","0").replace("Y","-1").replace("MT","-2"))

	if ch1_1==-1 or ch1_2==-1:
		line1=inputfile1.readline()
		continue

	this_category="somatic"

	start1_1=int(line1_split[4])
	start1_2=int(line1_split[7])

	stop1_1=int(line1_split[5])
	stop1_2=int(line1_split[8])

	orientation1=line1_split[9]

	while line2:
		line2_split=line2.rstrip().split("\t")
		ch2_1=int(line2_split[3].replace("X","0").replace("Y","-1").replace("MT","-2"))
		ch2_2=int(line2_split[6].replace("X","0").replace("Y","-1").replace("MT","-2"))
		if ch2_1==-1 or ch2_2==-1:
			line2=inputfile2.readline()
			continue
		orientation2=line2_split[9]

		if ch2_2 < ch1_2:
			line2=inputfile2.readline()
			continue
		if ch2_2 > ch1_2:
			cursor=inputfile2.tell()
			inputfile2.seek(-1*(min(cursor,300000)),1)
			line2=inputfile2.readline()
			line2=inputfile2.readline()
			break

		if orientation2!=orientation1:
			line2=inputfile2.readline()
			continue

		start2_1=int(line2_split[4])
		start2_2=int(line2_split[7])

		stop2_1=int(line2_split[5])
		stop2_2=int(line2_split[8])

		if stop2_2 < start1_2 - distance_thr:
			line2=inputfile2.readline()
			continue
		if start2_2 > stop1_2 + distance_thr:
			cursor=inputfile2.tell()
			inputfile2.seek(-1*(min(cursor,300000)),1)
			line2=inputfile2.readline()
			line2=inputfile2.readline()
			break
	
		#position 2 match
		
		if ch2_1 != ch1_1:
			line2=inputfile2.readline()
			continue
		if stop2_1 < start1_1 - distance_thr:
			line2=inputfile2.readline()
			continue
		if start2_1 > stop1_1 + distance_thr:
			line2=inputfile2.readline()
			continue
	
		#overlap
		this_category="match:"+line2_split[0]+":"+line2_split[1]+":"+line2_split[2]	
		cursor=inputfile2.tell()
		inputfile2.seek(-1*(min(cursor,300000)),1)
		line2=inputfile2.readline()
		line2=inputfile2.readline()
		break
	outputfile.write(line1.replace("\n","\t")+this_category+"\n")

	line1=inputfile1.readline()









