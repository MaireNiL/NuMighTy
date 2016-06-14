import sys
fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline()

ofn=fn.replace("r09_","r10_").replace(".ctr",".ctr.filtered.txt")
outputfile=file(ofn,"w")

toOutput=[]

while line:
	line_split=line.rstrip().split("\t")
	brkpoint=line_split[-1]

	if "breakpoints:" not in brkpoint:
		line=inputfile.readline()
		continue

	num_breakpoints=brkpoint.count(";")
	if num_breakpoints >6: #arbitrary
		line=inputfile.readline()
		continue
	brkpoint_split=brkpoint.replace("breakpoints:","").split(";")
	brkpoint1=brkpoint_split[0]
	base_overlap=int(brkpoint1.split(")-")[0].split("-(")[1])

	if base_overlap > 10:
		line=inputfile.readline()
		continue


	ch1_1=int(line_split[3].replace("X","0").replace("Y","-1").replace("MT","-2"))
	ch1_2=int(line_split[6].replace("X","0").replace("Y","-1").replace("MT","-2"))

	start1_1=int(line_split[4])
	stop1_1=int(line_split[5])

	start1_2=int(line_split[7])
	stop1_2=int(line_split[8])

	strand1_1=line_split[9].split(",")[0]
	strand1_2=line_split[9].split(",")[1]

	is_overlapped=0
	##CHECK any overlaps

	fn2=sys.argv[1]
	inputfile2=file(fn2)
	line2=inputfile2.readline()

	while line2:
		if line==line2: # same line
			line2=inputfile2.readline()
			continue
		line2_split=line2.rstrip().split("\t")
		brkpoint2=line2_split[-1]

		if "breakpoints:" not in brkpoint2:
			line2=inputfile2.readline()
			continue
		num_breakpoints2=brkpoint2.count(";")
		if num_breakpoints2 >6: # arbitrary
			line2=inputfile2.readline()
			continue

		ch2_1=int(line2_split[3].replace("X","0").replace("Y","-1").replace("MT","-2"))
		ch2_2=int(line2_split[6].replace("X","0").replace("Y","-1").replace("MT","-2"))

		start2_1=int(line2_split[4])
		stop2_1=int(line2_split[5])

		start2_2=int(line2_split[7])
		stop2_2=int(line2_split[8])

		strand2_1=line2_split[9].split(",")[0]
		strand2_2=line2_split[9].split(",")[1]

		if ch1_1!=ch2_2 or ch1_2!=ch2_1:
			line2=inputfile2.readline()
			continue
		if not(start1_1 >= start2_2 -200 and start1_1 <=start2_2 +200 and stop1_1 >=stop2_2 -200 and stop1_1 <=stop2_2 +200): #not overlap
			line2=inputfile2.readline()
			continue
		if not(start1_2 >= start2_1 -200 and start1_2 <=start2_1 +200 and stop1_2 >=stop2_1 -200 and stop1_2 <=stop2_1 +200): #not overlap
			line2=inputfile2.readline()
			continue
		if not(strand1_1 == strand2_2 and strand1_2 == strand2_1):
			line2=inputfile2.readline()
			continue
		#OVERLAP
		is_overlapped=1
		break

	pair_brk_info="not_detected"
	if is_overlapped==1:
		pair_brk_info=line2_split[1]+"::"+line2_split[-1]

	if ch1_1 >= ch1_2 :
		outputfile.write(line.replace("\n","\t")+pair_brk_info+"\n")
	else:
		if is_overlapped==0:
			outputfile.write(line.replace("\n","\t")+pair_brk_info+"\n")
	line=inputfile.readline()



print "done"
