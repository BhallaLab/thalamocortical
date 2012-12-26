#!/usr/bin/env python
# currentinjectiontest.py --- 

# 
# Filename: currentinjectiontest.py
# Description: 
# Author: 
# Maintainer: 
# Created: Fri Oct  5 10:25:32 2012 (+0530)
# Version: 
# Last-Updated: Wed Dec 26 09:38:48 2012 (+0530)
#           By: subha
#     Update #: 161
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# Run a current injection test on a specified cell
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

# from matplotlib import pyplot as plt
import numpy as np
import moose
from cell import TraubCell
from spinystellate import SpinyStellate
from suppyrRS import SupPyrRS
from suppyrFRB import SupPyrFRB
from supbasket import SupBasket
from supaxoaxonic import SupAxoaxonic
from supLTS import SupLTS
from tuftedIB import TuftedIB
from deepbasket import DeepBasket
from deepaxoaxonic import DeepAxoaxonic
from deepLTS import DeepLTS
from tuftedRS import TuftedRS
from nontuftedRS import NontuftedRS
from tcr import TCR
from nRT import nRT

simdt = 2.5e-6
plotdt = 1e-4

def setupmodel(celltype, delay, amplitude, duration):
    """Set up a single cell model with current injection of
    `amplitude` for `duration` period scheduled at `delay` time after
    start of simulation"""
    cellclass = eval(celltype)
    container = moose.Neutral('test_%s' % (celltype))
    model = moose.Neutral(container.path+'/model')
    cell = cellclass(cellclass.prototype, model.path+'/'+celltype)
    cell.method = 'hsolve'
    stim = moose.PulseGen(model.path+'/stim')
    stim.firstDelay = delay
    stim.firstWidth = duration
    stim.firstLevel = amplitude
    stim.secondDelay = 1e9
    stim.connect('outputSrc', cell.soma, 'injectMsg')
    data = moose.Neutral(container.path + '/data')
    vmtab = moose.Table(data.path+'/%s_soma_Vm' % (celltype))
    vmtab.stepMode = 3
    vmtab.connect('inputRequest', cell.soma, 'Vm')
    spiketab = moose.Table(data.path+'/%s_soma_spike' % (celltype))
    spiketab.stepMode = moose.TAB_SPIKE
    spiketab.stepSize = -20e-3
    spiketab.connect('inputRequest', cell.soma, 'Vm')
    return {'vm': vmtab,
	    'spike': spiketab,
	    'cell': cell,
	    'stim': stim,
	    'model': model,
	    'data': data}

def schedule(params):
    model = params['model']    
    data = params['data']
    moose.context.setClock(0, simdt)
    moose.context.setClock(1, simdt)
    moose.context.setClock(2, simdt)
    moose.context.setClock(3, simdt)
    moose.context.setClock(4, plotdt)
    moose.context.useClock(0, '%s/##[TYPE=Compartment]' % (model.path), 'init')    
    moose.context.useClock(1, '%s/##[TYPE=Compartment]' % (model.path), 'process')
    moose.context.useClock(2, '%s/##[TYPE=HSolve]' % (model.path), 'process')    
    moose.context.useClock(3, '%s/##[TYPE!=Compartment]' % (model.path), 'process')	
    moose.context.useClock(4, '%s/##' % (data.path), 'process')
    moose.context.reset()
    
def run_current_injection_test(celltype, delay, amplitude, duration, 
			       simtime, settlingtime=20e-3):
    components = setupmodel(celltype, delay, amplitude, duration)
    schedule(components)
    moose.context.step(simtime)
    vmdata = np.asarray(components['vm'])
    spikedata = np.asarray(components['spike'])
    tseries = np.linspace(0, simtime, len(vmdata))
    components['tseries'] = tseries
    return components

def save_data(components, fname):
    vmdata = np.transpose(np.vstack((components['tseries'], components['vm'])))
    vmfile = 'vm_'+fname
    np.savetxt(vmfile, vmdata)
    spikefile = 'spike_'+fname
    np.savetxt(spikefile, components['spike'])
    print 'Saved data in', vmfile, spikefile

import sys

default_delay = 100e-3
default_amplitude = 1e-9
default_duration = 100e-3
default_simtime = 300e-3
default_settlingtime = 20e-3

def print_usage(name):
	print 'Usage: %s celltype [delay amplitude duration simtime]' % (name)
	print 'Simulate for `simtime` a single instance of `celltype` \
with current injection of `amplitude` Ampere starting at `delay` \
second for `duration` time.'
	print 'If unspecified, delay=%g, amplitude=%g, duration=%g and \
simtime=%g second' % (default_delay, default_amplitude, 
		      default_duration, default_simtime)

if __name__ == '__main__':
    print sys.argv, len(sys.argv)
    if len(sys.argv) <= 1 or (len(sys.argv) > 2 and len(sys.argv) != 6):
	print_usage(sys.argv[0])
	sys.exit(2)
    celltype = sys.argv[1]
    delay = default_delay
    amplitude = default_amplitude
    duration = default_duration
    simtime = default_simtime
    settlingtime = default_settlingtime
    if len(sys.argv) == 6:
	delay = float(sys.argv[2])
	amplitude = float(sys.argv[3])
	duration = float(sys.argv[4])
	simtime = float(sys.argv[5])
    components = run_current_injection_test(
	celltype, delay, amplitude, duration, simtime, settlingtime)
    save_data(components, '%s_amp_%g.dat' % (celltype, amplitude))
    if len(np.nonzero(np.asarray(
		components['spike']) > settlingtime)[0]) > 0:
	sys.exit(0)
    sys.exit(1)


# 
# currentinjectiontest.py ends here
