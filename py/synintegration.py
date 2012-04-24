# synintegration.py --- 
# 
# Filename: synintegration.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Apr 24 15:30:05 2012 (+0530)
# Version: 
# Last-Updated: Tue Apr 24 21:25:41 2012 (+0530)
#           By: subha
#     Update #: 138
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This script is to set up a single spiny stellate cell and deliver
# artificial synaptic inputs to study the synaptic integration time. I
# noticed in network model that the spinystellate cells were
# responding about 40 ms after a stimulus. The purpose of this
# simulation is to find out the reason for that.
# 
# 

# Code:
import sys

sys.path.append('/data/subha/chamcham_moose/python')
import numpy as np
import h5py as h5
import moose
from trbsim import Simulation
from spinystellate import SpinyStellate

def get_syninfo(filename, cellname, srctype):
    """Load the synapses on cell from network file specified by
    filename"""
    ret = None
    with h5.File(filename, 'r') as netfile:
        syntab = netfile['/network/synapse'][:]
        indices = np.nonzero(np.char.startswith(syntab['dest'], cellname) & np.char.startswith(syntab['source'], srctype))[0]
        subarray = syntab[indices]
        print subarray
        comps = [int(entry[-1]) for entry in np.char.rpartition(subarray['dest'], '_')]
        print comps
        gbar = subarray['Gbar'].astype(float)
        tau1 = subarray['tau1'].astype(float)
        tau2 = subarray['tau2'].astype(float)
        Ek =  subarray['Ek'].astype(float)
        ret = (comps, subarray['type'], gbar, tau1, tau2, Ek)
    return ret

def test_synintegration(filename, target_cell, source_type):
    sim = Simulation('synintegration')
    cell = SpinyStellate(SpinyStellate.prototype, 
                         sim.model.path + '/SpinyStellate')
    vm_table = cell.comp[cell.presyn].insertRecorder('Vm', 'Vm', sim.data)
    # Create a common spike gen object
    sp = moose.SpikeGen('spike', sim.model)
    sp.threshold = 0.0
    sp.edgeTriggered = 1
    sptab = moose.Table('spike', sim.data)
    sptab.stepMode = 3
    sptab.connect('inputRequest', sp, 'state')
    (comp_indices, syntype, gbar, tau1, tau2, Ek, ) = \
        get_syninfo(filename, target_cell, source_type)    
    # Set up the synapses
    for ii in range(len(comp_indices)):
        print '%d\t%s\t%g\t%g\t%g\t%g' % (comp_indices[ii], 
                                          syntype[ii], 
                                          gbar[ii], 
                                          tau1[ii], 
                                          tau2[ii], 
                                          Ek[ii])
        comp = cell.comp[comp_indices[ii]]
        weight = 1.0
        if syntype[ii] == 'nmda':
            chan = moose.NMDAChan('nmda_from_%s' % (source_type), comp)
            chan.MgConc = 1.5            
            weight = gbar[ii]
        else:
            chan = moose.SynChan('%s_from_%s' % (syntype[ii], source_type), 
                                 comp)
        chan.Gbar = gbar[ii]
        chan.tau1 = tau1[ii]
        chan.tau2 = tau2[ii]
        chan.Ek = Ek[ii]
        comp.connect('channel', chan, 'channel')
        sp.connect('event', chan, 'synapse')
        count = chan.numSynapses
        if source_type == 'TCR':
            chan.delay[count-1] = 1e-3 # thalamocortical delay
        else:
            chan.delay[count-1] = 0.05e-3 # default delay
        chan.weight[count-1] = weight
        gktable = moose.Table('%s_%s_%s' % (cell.name, comp.name, chan.name), sim.data)
        gktable.stepMode = 3
        gktable.connect('inputRequest', chan, 'Gk')
    pulsegen = moose.PulseGen('pulse', sim.model)
    pulsegen.firstDelay = 3.0
    pulsegen.firstLevel = 1.0
    pulsegen.firstWidth = 1e-3
    pulsegen.trigMode = moose.FREE_RUN
    pulsegen.connect('outputSrc', sp, 'Vm')
    ptable = moose.Table('pulse', sim.data)
    ptable.stepMode = 3
    ptable.connect('inputRequest', pulsegen, 'output')
    sim.schedule(simdt=1e-6, plotdt=1e-6)
    sim.run(10.0)
    sim.save_data_h5(filename.replace('network_', 'synintegration_'))

if __name__ == '__main__':
    test_synintegration('/data/subha/rsync_ghevar_cortical_data_clone/2012_01_18/network_20120118_142820_7865.h5.new',
                        'SpinyStellate_0',
                        'TCR')
                        
            
# 
# synintegration.py ends here
