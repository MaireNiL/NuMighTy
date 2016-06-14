Example pipeline for extracing numts from bam00 //
Extract translocated reads with minimum mapping quality 10 from sam file generating bam00_tr_mq10.txt //
Cluster extracted reads using 6_clusterReadsOrientation.py generating bam00.clusters_ori.txt //
Run 7_filter.sh to generate bam00.clusters_ori_gt4.ctr //
Run 9_getSplitReads_local.py on r08_bam00_tr.somatic.ctr to generate r09_bam00_tr.somatic_brk1.ctr
Run 10_filterCluster.py on files r09_bam00_tr.somatic_brk1.ctr to generate r10_bam00_tr.somatic_brk1.ctr.filtered.txt 
