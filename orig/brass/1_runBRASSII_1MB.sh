#!/bin/bash
cp  `pipelineResult -db live -s samplenumber -p 814 -t 142`  samplenumber_brass2.tab ; python 01_getBRASSII_results_gt1mb.py samplenumber_brass2.tab
