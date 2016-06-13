import sys,math,random
fn=sys.argv[1]

inputfile=file(fn)
line=inputfile.readline()

ofn=fn.replace(".txt",".circ")
outputfile=file(ofn,"w")
outputfile.write("# chr	start	end	value\n")

prev_ch=""

available_chr=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT","M"]
while line:
	line_split=line.rstrip().split("\t")

	ch=line_split[6]
	if ch not in available_chr:
		line=inputfile.readline()
		continue

	pos=int(line_split[7])

	if ch!=prev_ch:
		prev_ch=ch
		chr_indice=0
	chr_indice+=1
	if chr_indice==1:
		pos_diff=10000000
	else:
		pos_diff=pos-prev_pos
		if pos_diff==0:
			pos_diff=random.uniform(0.1,1.0)
	pos_diff_log=math.log(pos_diff,10)
	outputfile.write("hs"+ch+"\t"+str(pos)+"\t"+str(pos+100)+"\t"+str(round(pos_diff_log,5))+"\n")
	prev_pos=pos

	line=inputfile.readline()


