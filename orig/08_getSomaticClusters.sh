#!/bin/bash
# compare clusters from tumour sample against paired normals
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4.ctr 01_paired_normals/00006_read_clusters/hostsample00_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn.ctr 01_paired_normals/00006_read_clusters/hostsample01_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample02_tr.clusters_ori_s.ctr
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample03_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample04_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample05_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample06_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample07_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample08_tr.clusters_ori_s.ctr 
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample09_tr.clusters_ori_s.ctr
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample10_tr.clusters_ori_s.ctr
grep -cv "match" sample00*
python 08_compare_paired_normal.py sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn.ctr 01_paired_normals/00006_read_clusters/hostsample11_tr.clusters_ori_s.ctr
grep -cv "match" sample00*
grep -v "match" sample00_tr_mq10.clusters_ori_gt4_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn_pn.ctr > sample00_tr_mq10.somatic.ctr
