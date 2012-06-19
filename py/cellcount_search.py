# cellcount_search.py --- 
# 
# Filename: cellcount_search.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jun 18 14:29:18 2012 (+0530)
# Version: 
# Last-Updated: Mon Jun 18 15:42:08 2012 (+0530)
#           By: subha
#     Update #: 107
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This script runs multiple simulations by updating the custom.ini
# file with different cell counts.
# 
# 

# Change log:
# 
# 
# 
# 

# Code:

import ConfigParser as cp
import os
import random
import uuid
import subprocess as sp
from datetime import datetime
from collections import deque

starting_cellcount = {
    'SupPyrRS': 100,
    'SupPyrFRB': 5,
    'SupBasket': 9 ,
    'SupAxoaxonic': 9,
    'SupLTS': 9,
    'SpinyStellate':	240,
    'TuftedIB': 0,
    'TuftedRS': 0,
    'DeepBasket': 0,
    'DeepAxoaxonic': 0,
    'DeepLTS': 0,
    'NontuftedRS': 0,
    'TCR': 100,
    'nRT': 0 
    }

# Keep increasing the candidate cell populations by `increments` for `steps` steps.
candidates = ['DeepBasket', 'DeepAxoaxonic', 'DeepLTS']
increments = 5
steps = 11
process_count = 5

if __name__ == '__main__':
    # Store the original configuration ahead of changing the file.
    config = cp.SafeConfigParser()
    config.read(['custom.ini'])
    # tmp_config is what we'll keep changing
    tmp_config = cp.SafeConfigParser()
    tmp_config.read(['custom.ini'])    
    for celltype, count in starting_cellcount.items():
        tmp_config.set('cellcount', celltype, str(count))
    processes = deque()
    files = deque()
    for celltype in candidates:
        cellcount = int(tmp_config.get('cellcount', celltype))
        notes = 'Cell count of %s: ' % (celltype)
        for step in range(steps):
            cellcount += step * increments
            notes = '%s %d.' % (notes, cellcount)
            tmp_config.set('cellcount', celltype, str(cellcount)) 
            net_rng = random.randint(0, 0xffffffff)
            tmp_config.set('numeric', 'numpy_rngseed', str(net_rng))
            with open('custom.ini', 'w') as configfile:
                tmp_config.write(configfile)
            notes = '%s numpy RNG seed: %d' % (notes, net_rng)
            for pcount in range(process_count):                
                timestamp = datetime.now()
                simid = uuid.uuid4().int
                outfilename = 'trbsim_%s_%d.out' % (timestamp.strftime('%Y%m%d_%H%M%S'), simid)
                errfilename = 'trbsim_%s_%d.err' % (timestamp.strftime('%Y%m%d_%H%M%S'), simid)                
                outfile = open(outfilename, 'w')
                errfile = open(errfilename, 'w')
                proc = sp.Popen(['python', 'trbsim.py', '-n', notes], stdout=outfile, stderr=errfile)
                print 'Starting PID', proc.pid, 'ouput in', outfilename
                processes.append(proc)
                files.append((outfile, errfile))
            while len(processes) > 0:
                p = processes.popleft()
                p.wait()
                f = files.popleft()
                f[0].close()
                f[1].close()


# 
# cellcount_search.py ends here
