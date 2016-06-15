echo sample000 ; samtools view -F 8 -q 10  sample000.bam | awk '$7!="="' | sort -k7,7n -k8,8n > sample000_tr_mq10.txt  &
echo sample001 ; samtools view -F 8 -q 10  sample001.bam | awk '$7!="="' | sort -k7,7n -k8,8n > sample001_tr_mq10.txt  &
echo sample002 ; samtools view -F 8 -q 10  sample002.bam | awk '$7!="="' | sort -k7,7n -k8,8n > sample002_tr_mq10.txt  &
