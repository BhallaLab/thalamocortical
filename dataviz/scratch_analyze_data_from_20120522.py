# scratch.py --- 
# 
# Filename: scratch.py
# Description: 
# Author: 
# Maintainer: 
# Created: Fri Jun 29 09:13:54 2012 (+0530)
# Version: 
# Last-Updated: Fri Jul  6 09:08:03 2012 (+0530)
#           By: subha
#     Update #: 131
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
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

if __name__ == '__main__':
    df = h5.File('/data/subha/cortical/py/data/2012_05_22/data_20120522_152734_10973.h5', 'r')
    schedict = dict(df['runconfig/scheduling'])
    simtime = float(schedict['simtime'])
    stimdict = dict(df['runconfig/stimulus'])
    onset = float(stimdict['onset']) + float(stimdict['bg_interval'])
    print 'Stimulus onset:', onset
    bg = df['/stimulus/stim_bg'][:]
    dt = simtime / len(bg)
    bgtimes = np.nonzero(np.diff(bg) > 0)[0] * dt
    bgtimes = np.r_[bgtimes, simtime]
    print 'Background times', bgtimes
    bins = zip(bgtimes[:-1], bgtimes[1:])
    bgbins = [bins[ii] for ii in range(0, len(bins), 2)]
    print 'Background bins'
    print bgbins    
    probebins = [bins[ii] for ii in range(1, len(bins), 2)]
    print 'Probe bins',
    print probebins
    spikedict = dict([(cell, df['/spikes/%s' % (cell)][:]) for cell in df['/spikes'] if cell.startswith('Spiny')])
    binnedbg = defaultdict(list)
    binnedprobe = defaultdict(list)
    bad_count = defaultdict(list)
    for cell, spiketimes in spikedict.items():
        binnedbg[cell] = [spiketimes[(spiketimes >= bin[0]) & (spiketimes < bin[1])] - bin[0] for bin in bgbins]
        binnedprobe[cell] = [spiketimes[(spiketimes >= bin[0]) & (spiketimes < bin[1])] - bin[0] for bin in probebins]
    fig = plt.figure()
    ax_spikes = fig.add_subplot(2, 1, 1, title='%s\nSpike times wrt stimulus\nBlue x: background only, Red +: background + probe' % df.name)
    ax_spikes.set_xlabel('Time from stimulus')
    ax_spikes.set_ylabel('Cell #')
    ax_psth = fig.add_subplot(2, 1, 2, title='PSTH')
    ax_psth.set_xlabel('Spike time wrt stimulus onset (2 ms time bin)')
    ax_psth.set_ylabel('Spike count')
    timebins = np.arange(0, bgbins[0][1] - bgbins[0][0], 0.002)
    bghist = np.zeros(len(timebins)-1)
    probehist = np.zeros(len(timebins)-1)    
    # ax_counts = fig.add_subplot(2, 1, 2, title='Spike count follwoing stimulus\nBlue x: background only, Red +: background + probe')    
    # ax_counts.set_xlabel('Stimulus #')
    # ax_counts.set_ylabel('Number of spikes from a cell')
    for cell, spiketimes in binnedbg.items():
        cellno = int(cell.rpartition('_')[-1])
        spikes = np.concatenate(spiketimes)
        ax_spikes.plot(spikes, np.ones(len(spikes)) * cellno, 'b,')
        hist, edges = np.histogram(spikes, bins=timebins)
        bghist = bghist + hist
        # counts = np.array([len(bin) for bin in spiketimes])
        # ax_counts.plot(np.arange(0, len(counts), 1.0), counts, 'bx')
    for cell, spiketimes in binnedprobe.items():
        cellno = int(cell.rpartition('_')[-1])
        spikes = np.concatenate(spiketimes)
        ax_spikes.plot(spikes, np.ones(len(spikes)) * cellno, 'r,')
        hist, edges = np.histogram(spikes, bins=timebins)
        probehist += hist
        # counts = np.array([len(bin) for bin in spiketimes])
        # ax_counts.plot(np.arange(0, len(counts), 1.0), counts, 'r+')
    ax_psth.plot((timebins[:-1] + timebins[1:])/2.0, bghist, 'bv-.')
    ax_psth.plot((timebins[:-1] + timebins[1:])/2.0, probehist, 'r^-.')
    plt.show()

# 
# scratch.py ends here
