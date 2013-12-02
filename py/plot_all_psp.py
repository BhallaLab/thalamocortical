# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Tue Jan 22 11:51:15 2013 (+0530)
#           By: subha
#     Update #: 457
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
# Thu Dec 27 16:54:13 IST 2012 - adding function to test a single
# spiny stellate with TCR cells converging onto it.
# 
# 

# Code:

import numpy as np
import pylab
import matplotlib
from matplotlib import pyplot as plt
import sys
sys.path.append('.')
# sys.path.append('/data/subha/chamcham_moose/python')

import moose
print moose.__file__
from moose import _moose
print _moose.__file__

from cell import *
from simulation import Simulation
from trbnetdata import TraubFullNetData

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

import synapse
import random
from itertools import cycle, izip, chain
from compartment import compare_compartment
import config
import os
from datetime import datetime

netdata = TraubFullNetData()

config.solver = 'hsolve'

def simulate_psp(simdt, plotdt, simtime):
    sim = Simulation('dump_psp')
    for pretype in netdata.celltype:
        for posttype in netdata.celltype:
            for chantype in ['ampa', 'nmda', 'gaba']:                
                setup(pretype, posttype, chantype, sim)
    sim.simdt = simdt
    sim.plotdt = plotdt
    print 'Scheduling'
    sim.schedule()
    print 'Starting simulation'
    sim.run(simtime)
    print 'Simulation over'
    datadir = 'all_pair_psp'
    try:
        os.mkdir(datadir)
    except OSError, e:
        print e
    for tab in config.context.getWildcardList('%s/##[TYPE=Table]' % (sim.data.path), True):
        ts = np.linspace(0, simtime, len(tab))
        data = np.vstack((ts, data)).transpose()
        dpath = os.path.join(datadir, tab.parent.name)
        try:
            os.mkdir(dpath)
        except OSError, e:
            print e                    
        fname = '%s/%s.dat' % (dpath, tab.name)
        np.savetxt(fname, data)
        if 'postVm' in tab.name:
            print '##### ',  tab.path, ': max=', np.max(tab)
        print 'Saved data from %s to %s' % (tab.name, fname)
    print 'Finished saving data'
    
def setup(pretype, posttype, chantype, sim):
    """Plot the PSP in`posttype` cell from `chantype` synapse due to spiking in `pretype` cell.

    `chantype` can be `nmda`, `ampa` and `gaba`"""    
    preclass = eval(pretype)
    postclass = eval(posttype)
    preidx = netdata.celltype.index(pretype)
    postidx = netdata.celltype.index(posttype)
    # We arbitrarily choose the first entry in allowed compartment
    # list - since the list is ordered, this is likely to be the one
    # closest to soma.
    try:
        postcomp = netdata.allowed_comps[preidx][postidx][0]
    except IndexError:
        print 'No synapses from %s to %s' % (pretype, posttype)
        return
    Ek = 0
    tau1 = 0
    tau2 = 1.0
    gbar = 1.0
    if chantype == 'ampa':
        gbar = netdata.g_ampa_baseline[preidx][postidx]
        tau1 = netdata.tau_ampa[preidx][postidx]
        tau2 = tau1
    elif chantype == 'nmda':
        gbar = netdata.g_nmda_baseline[preidx][postidx]
        tau1 = netdata.tau_nmda[preidx][postidx]
        tau2 = 5e-3
    else:
        Ek = netdata.ek_gaba[postidx]
        gbar = netdata.g_gaba_baseline[preidx][postidx]
        tau1 = netdata.tau_gaba[preidx][postidx]
        tau2 = 0.0
    if gbar <= 0.0:
        print 'No %s synapse from %s to %s' % (chantype, pretype, posttype)
        return 
    model_container = moose.Neutral('%s/%s_%s_%s' % (sim.model.path, chantype, pretype, posttype))
    precell = preclass(preclass.prototype, '%s/%s_pre' % (model_container.path, pretype))
    postcell = postclass(postclass.prototype, '%s/%s_post' % (model_container.path, posttype))
    syn = precell.comp[precell.presyn].makeSynapse(postcell.comp[postcomp],
                                                   name='%s_%s_%s' % (chantype, pretype, posttype),
                                                   Ek=Ek,
                                                   Gbar=gbar,
                                                   tau1=tau1, 
                                                   tau2=tau2,
                                                   delay=1e-3)
    if chantype == 'nmda':
        syn.MgConc = netdata.MgConc
        syn.saturation = 1.0
    stim = moose.PulseGen('%s/stimulus' % (model_container.path))
    stim.delay[0] = 2.0
    stim.level[0] = 1e-9
    stim.width[0] = 10e-3
    stim.connect('outputSrc', precell.soma, 'injectMsg')
    data_container = moose.Neutral('%s/%s_%s_%s' % (sim.data.path, chantype, pretype, posttype))
    datatables = []
    datatables.append(moose.Table('%s/preVm_%s_%s_%s' % (data_container.path, chantype, pretype, posttype)))
    datatables[-1].connect('inputRequest', precell.soma, 'Vm')
    datatables.append(moose.Table('%s/postVm_%s_%s_%s' % (data_container.path, chantype, pretype, posttype)))
    datatables[-1].connect('inputRequest', postcell.soma, 'Vm')
    datatables.append(moose.Table('%s/synIk_%s_%s_%s' % (data_container.path, chantype, pretype, posttype)))
    datatables[-1].connect('inputRequest', syn, 'Ik')
    for tab in datatables:
        tab.stepMode = 3

if __name__ == '__main__':
    # test_tcr_spinstell_ampa()
    # test_tcr_ss_spiking(0.0)
    simulate_psp(5e-6, 1e-4, 3.0)
# 
# test_tcr_spinstell.py ends here