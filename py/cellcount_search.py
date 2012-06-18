# cellcount_search.py --- 
# 
# Filename: cellcount_search.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jun 18 14:29:18 2012 (+0530)
# Version: 
# Last-Updated: Mon Jun 18 15:22:32 2012 (+0530)
#           By: subha
#     Update #: 83
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
import subprocess as sp
from datetime import datetime
from collections import deque

startiong_cellcount = {
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
processes = 5

if __name__ == '__main__':
    # Store the original configuration ahead of changing the file.
    config = cp.SafeConfigParser()
    config.read(['custom.ini'])
    # tmp_config is what we'll keep changing
    tmp_config = cp.SafeConfigParser()
    tmp_config.read(['custom.ini'])    
    processes = deque()
    files = deque()
    for celltype in candidates:
        cellcount = int(tmp_config.get('cellcount', celltype))
        notes = 'Cell count of %s: ' % (celltype)
        for step in range(steps):
            cellcount += step * increments
            notes = '%s %d.' % (notes, cellcount)
            tmp_config.set('cellcount', celltype, str(cellcount)) 
            net_rng = random.random(0, 0xffffffff)
            tmp_config.set('numeric', 'numpy_rngseed', str(net_rng))
            tmp_config.write('custom.ini')
            notes = '%s numpy RNG seed: %d' % (notes, net_rng)
            for pcount in range(processes):                
                timestamp = datetime.now()
                outfilename = 'trbsim_%s.out' % (timestamp.strftime('%Y%m%d_%H%M%s'))
                errfilename = 'trbsim_%s.err' % (timestamp.strftime('%Y%m%d_%H%M%s'))                
                outfile = open(outfilename, 'w')
                errfile = open(errfilename, 'w')
                proc = sp.Popen(['python', 'trbsim.py', '-n', notes], stdout=outfile, stderr=errfile)
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
