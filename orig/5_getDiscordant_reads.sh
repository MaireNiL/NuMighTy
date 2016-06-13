#!/bin/bash
echo sample00_MT.bam; samtools view -F 14 sample00_MT.bam | awk '$7=="=" && $2!=145 && $2!=81 && $2!=161 && $2!=97' | wc -l
echo sample01_MT.bam; samtools view -F 14 sample01_MT.bam | awk '$7=="=" && $2!=145 && $2!=81 && $2!=161 && $2!=97' | wc -l
