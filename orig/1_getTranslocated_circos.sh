#!/bin/bash
echo 0 ; samtools view -F 8 sample0_MT.bam | awk '$7!="=" && $7!="Y"' | sort -k7,7n -k8,8n > sample0_MTtr.txt ; python 2_forCircos.py sample0_MTtr.txt
echo 1 ; samtools view -F 8 sample1.bam | awk '$7!="=" && $7!="Y"' | sort -k7,7n -k8,8n > sample1_MTtr.txt ; python 2_forCircos.py sample1_MTtr.txt
