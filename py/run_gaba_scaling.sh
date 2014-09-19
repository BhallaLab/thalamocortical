#!/bin/bash
moosedir=/home/subhasis/moose_trunk_20131016/python
pypath="$moosedir:$PYTHONPATH"
cd /home/subhasis/cortical/py
# Line # 77 contains GABA conductance_scale
# sed -i '77s/conductance_scale.*$/conductance_scale = 1.5/'  custom.ini
# for hostindex in `seq 40 49`; do
# 	hostid=compute-0-"$hostindex"
#       echo "Starting job on $hostid"
#       nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath trbsim.py -n "Running with GABA-scale=1.5" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
#       sleep 3
sed -i '77s/conductance_scale.*$/conductance_scale = 1.2/'  custom.ini
for hostindex in `seq 30 39`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=1.2" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done      
sed -i '77s/conductance_scale.*$/conductance_scale = 1.1/'  custom.ini
for hostindex in `seq 30 39`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=1.1" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done      
sed -i '77s/conductance_scale.*$/conductance_scale = 0.9/'  custom.ini
for hostindex in `seq 20 29`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=0.9" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done      
sed -i '77s/conductance_scale.*$/conductance_scale = 0.8/'  custom.ini
for hostindex in `seq 0 9`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=0.8" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done      
sed -i '77s/conductance_scale.*$/conductance_scale = 0.7/'  custom.ini
for hostindex in `seq 40 49`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=0.7" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done
sed -i '77s/conductance_scale.*$/conductance_scale = 0.6/'  custom.ini
for hostindex in `seq 40 49`; do
  hostid=compute-0-"$hostindex"
      echo "Starting job on $hostid"
      nohup mpirun -np 1 --host $hostid -x LD_LIBRARY_PATH -x PYTHONPATH=$pypath /usr/local/bin/python2.7 trbsim.py -n "Running with GABA-scale=0.6" &> nohup_"$hostid"_`date '+%Y%m%d_%H%M%S'`.log &
      sleep 3
done

