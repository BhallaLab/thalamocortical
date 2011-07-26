#!/bin/bash
#nohup python trbsim.py -t 1.0 --reseed | gzip -c - > nohup_`date '+%Y%m%d_%H%M%S'`.gz 
#python trbsim.py -t 1.0 --reseed 2>&1 | gzip -c - > nohup_`date '+%Y%m%d_%H%M%S'`.gz
at now + 36 hours -f /data/subha/cortical/py/commandline.sh
pushd /data/subha/cortical/py/
python26 trbsim.py -t 5.0 -x 0.8 --reseed --stochastic > nohup_`date '+%Y%m%d_%H%M%S'` 2>&1 &
popd
