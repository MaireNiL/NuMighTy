#!/bin/bash
cp  `pipelineResult -db live -s samplenumber -p 814 -t 142`  samplenumber_brass2.tab ; python 1_getBRASSII_results.py samplenumber_brass2.tab
