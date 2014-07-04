#!/bin/bash
count=10
if (($# > 0)); then
    count=$1
fi
dbcount=40
if (($# > 1)); then
    dbcount=$2
fi
distr=unknown
if (($# > 2)); then
    distr=$3
fi

for ((ii=1; ii <= $count; ++ii)); do
    echo "$ii"
    nohup python trbsim.py -n "Running #$ii-th simulation of $dbcount deep basket cells with $distr distribution of synaptic conductances." &> nohup_`date '+%Y%m%d_%H%M%S'`.log &
    wait ${!}
    echo "Done running #$ii-th simulation with $dbcount deep basket cells"
done


