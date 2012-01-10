# search_tcr_isi.py --- 
# 
# Filename: search_tcr_isi.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Jan 10 17:21:25 2012 (+0530)
# Version: 
# Last-Updated: Tue Jan 10 18:48:23 2012 (+0530)
#           By: subha
#     Update #: 57
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# 
# 
# 

# Change log:
# 
# 
# 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 
# 

# Code:

import os
import subprocess
from datetime import datetime
import numpy
import h5py
import sys
start = 5e-3
stop = 500e-3
dt = 5e-3
if __name__ ==  '__main__':
    data_dir = os.path.join('data', datetime.now().strftime('%Y_%m_%d'))
    last_isi = 1e10
    current_isi = start
    while current_isi < stop:
        proc = subprocess.Popen(['python2.6', 'trbsim.py', '-i', str(current_isi), '-n', 'running with isi=%g' % (current_isi)])
        pid = proc.pid
        print 'PID', pid
        proc.wait()
        for filename in os.listdir(data_dir):
            if filename.find(str(pid)) > 0 and filename.startswith('data'):
                print 'Checking', filename
                datafile = h5py.File(os.path.join(data_dir, filename), 'r')
                spike_times = numpy.array(datafile['/spikes/TCR_0'])
                if len(spike_times) > 2 and (spike_times[-1] - spike_times[-2]) < 1.5 * current_isi:
                    print 'ISI:', current_isi, 'works'
                    print spike_times
                    datafile.close()
                    sys.exit(0)
                print 'ISI:', current_isi, 'does not work'
                print spike_times
        current_isi += dt
                    
            


# 
# search_tcr_isi.py ends here
