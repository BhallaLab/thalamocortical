#!/bin/bash
#nohup python trbsim.py -t 1.0 --reseed | gzip -c - > nohup_`date '+%Y%m%d_%H%M%S'`.gz 
#python trbsim.py -t 1.0 --reseed 2>&1 | gzip -c - > nohup_`date '+%Y%m%d_%H%M%S'`.gz
at now + 36 hours -f /data/subha/cortical/py/commandline.sh
pushd /data/subha/cortical/py
outfile=nohup_`date '+%Y%m%d_%H%M%S'`

#export PYTHONPATH="$PYTHONPATH":"/usr/local/lib/python2.6/site-packages":"/usr/local/lib/python2.6/site-packages/numpy/core/"
echo "PYTHONPATH=$PYTHONPATH" >> $outfile
echo "SHELL=$SHELL" >> $outfile
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/usr/local/lib
#python2.6 trbsim.py -c cellcount.csv -t 0.005 -x 0.8 --reseed --stochastic >> $outfile 2>&1 &
python2.6 trbsim.py -d 0.25e-3 -p 1e-3 -t 5.0 -x 0.8 --reseed --stochastic >> $outfile 2>&1 &
popd
