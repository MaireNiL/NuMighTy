Example pipeline for extracing numts from bam_00 \\
Extract translocated reads with minimum mapping quality 10 from sam file generating bam_00_tr_mq10.txt \\
Cluster extracted reads using 6_clusterReadsOrientation.py generating bam_00.clusters_ori.txt \\
Run 7_filter.sh to generate bam_00.clusters_ori_gt4.ctr \\
