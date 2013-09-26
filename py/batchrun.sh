#!/bin/bash
count=10
if (($# > 0)); then
    count=$1
fi
dbcount=55
if (($# > 1)); then
    dbcount=$2
fi
for ((ii=1; ii <= $count; ++ii)); do
    echo "$ii"
    nohup python trbsim.py -n "Running #$ii-th simulation of $dbcount deep basket cells with normal distribution of synaptic conductances but no variability in individual cells." &> nohup_`date '+%Y%m%d_%H%M%S'`.log &
    wait ${!}
    echo "Done running #$ii-th simulation with $dbcount deep basket cells"
done


