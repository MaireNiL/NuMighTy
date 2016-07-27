Example pipeline for extracing numts from bam00 
## Extract translocated reads with minimum mapping quality 10 from sam file 
bash 1_getTranslocated_reads.sh bam00.bam > bam00_tr_mq10.txt 
## Collapse translocated reads into clusters
# 6_clusterReadsOrientation.py allows cluster with single supporting read
python 6_clusterReadsOrientation.py r06_bam00_tr_mq10.txt > r07_bam00_tr.clusters_ori.ctr
## Sort paired normal clusters 
bash 7_sortPairedNormal.sh r07_bam00_tr.clusters_ori.ctr > r07_bam00_tr.clusters_ori_s.ctr 
## Find somatic clusters with more than 4 supporting reads
bash 7_filter.sh bam00.clusters_ori_gt4.ctr 
or just 
bash 7_filter.sh #extracts clusters supported by greater than for reads to generate *clusters_ori_gt4.ctr file from *clusters_ori.ctr
bash 8_getSomaticClusters.sh #runs 8_compare_paired_normal.py to generate *.somatic.ctr from *.clusters_ori_gt4.ctr
python 9_getSplitReads_local.py r08_bam00_tr.somatic.ctr > r09_bam00_tr.somatic_brk1.ctr
python 10_filterCluster.py r09_bam00_tr.somatic_brk1.ctr > r10_bam00_tr.somatic_brk1.ctr.filtered.txt 
