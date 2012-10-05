# currentinjectiontestdriver.py --- 
# 
# Filename: currentinjectiontestdriver.py
# Description: 
# Author: 
# Maintainer: 
# Created: Fri Oct  5 13:58:15 2012 (+0530)
# Version: 
# Last-Updated: Fri Oct  5 17:38:47 2012 (+0530)
#           By: subha
#     Update #: 48
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

import sys
import subprocess
import numpy as np

default_delay = 100e-3
default_duration = 100e-3
default_simtime = 300e-3

if __name__ == '__main__':
    celltype = sys.argv[1]
    logamps= np.linspace(-11, -9, 5)
    amps = np.power(10.0, logamps)
    for amp in amps:
	proc = subprocess.Popen('./currentinjectiontest.py %s %g %g %g %g' %
				(celltype, 
				 default_delay,
				 amp, 
				 default_duration, 
				 default_simtime),
				stdin=subprocess.PIPE,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE, shell=True)
	out, err = proc.communicate()
	print '*Output for amp = %g*' % (amp)
	print out
	print '*Error for amp = %g*' % (amp)
	print err
	print '---'
	ret = proc.poll()
	if ret == 0:
	    print 'spiking at %g A' % (amp)
	    sys.exit(0)
    print 'No spiking in range: %g - %g' % (amps[0], amps[-1])


# 
# currentinjectiontestdriver.py ends here
