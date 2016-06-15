#!/bin/bash
echo bam00_MT.bam; samtools view -F 14 bam00_MT.bam | awk '$7=="=" && $2!=145 && $2!=81 && $2!=161 && $2!=97' | wc -l
echo bam01_MT.bam; samtools view -F 14 bam01_MT.bam | awk '$7=="=" && $2!=145 && $2!=81 && $2!=161 && $2!=97' | wc -l
