# analyzer_bursts.py --- 
# 
# Filename: analyzer_bursts.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Aug  6 09:10:30 2012 (+0530)
# Version: 
# Last-Updated: Tue Aug  7 09:31:11 2012 (+0530)
#           By: subha
#     Update #: 105
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

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    strdata = np.loadtxt('filtered.txt', dtype=str)
    # dtype will be :
    #
    # [('simtime', float), ('bginterval', float), ('ppinterval',
    # float), ('spikecount', float), ('cellcount', float),
    # ('inhibitory', float), ('tcr', float), ('stimulated', float),
    # ('burstlength', float), ('spikesperburst', float)]
    dtype = [(col, float) for col in strdata[0, 1:]]
    numdata = np.array(strdata[1:,1:], dtype=float)
    numdata = numdata.view(dtype).copy()
    idx_bursting = np.nonzero(numdata['spikesperburst'])[0]
    print 'file inhibitory# spikes_per_burst'
    for ii in idx_bursting:
        print strdata[ii+1, 0], strdata[ii+1, 6], strdata[ii+1, -1]
    idx_nonbursting = np.nonzero(numdata['spikesperburst'] == 0)[0]
    plt.subplot(211)
    plt.ylabel('no. of simulations')
    plt.xlabel('inhibitory cell count')
    counts, bins, patches = plt.hist([numdata[idx_bursting]['inhibitory'], numdata[idx_nonbursting]['inhibitory']], label=['Bursting', 'Not-bursting'])
    plt.xticks(bins)
    plt.legend()
    plt.subplot(212)
    plt.xlabel('inhibitory cell count / spiny stellate cell count')
    plt.ylabel('no. of simulations')
    exc_inh = numdata['inhibitory'] / numdata['cellcount']
    outlier = np.nonzero((exc_inh > 1.0) & (numdata['spikesperburst'] > 0))[0]
    if outlier < len(strdata):
        print 'OUTLIER - bursting in spite of high inhibition'
        print strdata[outlier+1]
    # Now let us check which simulations are showing no bursts in
    # spite of exc/inh ratio being greater than that of some bursting
    # cells.
    bursting = np.nonzero(numdata['spikesperburst'] > 0)[0]
    nonbursting = np.nonzero(numdata['spikesperburst'] <= 0)[0]
    bdata = []
    nbdata = []
    if len(bursting) > 0:
        bfiles = strdata[bursting+1, 0]
        bdata = exc_inh[bursting].copy()
    if len(nonbursting) > 0:
        nbfiles = strdata[nonbursting+1, 0]
        nbdata = exc_inh[nonbursting].copy()
    if (len(bursting) > 0) and (len(nonbursting) > 0):
        outliers = np.nonzero(nbdata <= max(bdata))[0]
        print 'Nonbursting with too little inhibition:'
        for index in outliers:
            idx = np.where(strdata[:,0] == nbfiles[index])
            for ii in idx:
                print strdata[ii]
    
    counts, bins, patches = plt.hist([exc_inh[idx_bursting], exc_inh[idx_nonbursting]], label=['bursting', 'non-bursting'])
    plt.xticks(bins)
    plt.legend()
    plt.show()
    


# 
# analyzer_bursts.py ends here
