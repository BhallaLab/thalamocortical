# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Tue Jan 22 11:49:05 2013 (+0530)
#           By: subha
#     Update #: 456
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
sys.path.append('/data/subha/chamcham_moose/python')

import moose
print moose.__file__
from moose import _moose
print _moose.__file__
from simulation import Simulation
from trbnetdata import TraubFullNetData
from spinystellate import SpinyStellate
from tcr import TCR
import synapse
import random
from itertools import cycle, izip, chain
from compartment import compare_compartment
import config
import os
from datetime import datetime

def test_tcr_ss_spiking(offset=0.0):
    """This for checking the effect of TCR cells on SpinyStellate (SS)
    cells.

    I want to see what is the minimum number of presynaptic TCR cells
    needed to spike within a short interval to stimulate a
    post-synaptic SS. Depending on the value of {offset}, the TCR
    cells may receive stimulus at slightly different times (differing
    by offset).

    I create an array of SS cells and an array of TCR cells. SS_n
    receives input from TCR_0 trough TCR_n.

    """
    print 'test_tcr_ss_spiking: starting at', datetime.now().strftime('%Y%m%d_%H%M%S')
    datadir = 'data_%s' % (datetime.now().strftime('%Y%m%d_%H%M%S'))
    os.mkdir(datadir)
    netdata = TraubFullNetData()
    config.solver = 'hsolve'
    sim = Simulation('tcr_ss')
    tcr_idx = netdata.celltype.index('TCR')
    ss_idx = netdata.celltype.index('SpinyStellate')
    num_ss, num_tcr = 19, 20
    ss = [SpinyStellate(SpinyStellate.prototype, '%s/SS_%d' % (sim.model.path, idx)) for idx in range(num_ss)]
    tcr = [TCR(TCR.prototype, '%s/TCR_%d' % (sim.model.path, idx)) for idx in range(num_tcr)]
    print 'test_tcr_ss_spiking: cells created at', datetime.now().strftime('%Y%m%d_%H%M%S')
    pre_per_post = netdata.pre_post_ratio[tcr_idx][ss_idx]
    nmda_tabs = []
    ampa_tabs = []
    vm_tabs = []
    ca_tabs = []
    for index, cell in enumerate(ss):
        # Select the compartments for creating synapses. On n-th cell
        # we create n+1 synapses (counting from 0).
        post_comp_list = [cell.comp[ii] for ii in random.sample(netdata.allowed_comps[tcr_idx][ss_idx], index+1)]
        print cell.path, 'receiving input on %d comps:' % (index+1)
        # The first n+1 TCR cells will send input to n-th spiny stellate cell.
        for precell, postcomp in zip(tcr, post_comp_list):
            print '\t', precell.path, 'on', postcomp.name
            ampa = precell.comp[precell.presyn].makeSynapse(postcomp,
                                                name='ampa__%s__%s__%s' % (precell.name, cell.name, postcomp.name),
                                                classname='SynChan',
                                                Ek=0.0,
                                                Gbar=netdata.g_ampa_baseline[tcr_idx][ss_idx],
                                                tau1=netdata.tau_ampa[tcr_idx][ss_idx],
                                                tau2=netdata.tau_ampa[tcr_idx][ss_idx],
                                                delay=1e-3)
            ampa_tabs.append(moose.Table('%s/g%s' % (sim.data.path, ampa.name)))
            ampa_tabs[-1].stepMode = 3
            ampa_tabs[-1].connect('inputRequest', ampa, 'Gk')
            nmda = precell.comp[precell.presyn].makeSynapse(postcomp,
                                                name='nmda__%s__%s__%s' % (precell.name, cell.name, postcomp.name),
                                                classname='NMDAChan',
                                                Ek=0.0,
                                                Gbar=netdata.g_nmda_baseline[tcr_idx][ss_idx],
                                                tau1=netdata.tau_nmda[tcr_idx][ss_idx],
                                                tau2=5e-3,
                                                delay=1e-3)
            nmda.MgConc = netdata.MgConc
            nmda.saturation = 1.0
            nmda_tabs.append(moose.Table('%s/g%s' % (sim.data.path, nmda.name)))
            nmda_tabs[-1].stepMode = 3
            nmda_tabs[-1].connect('inputRequest', nmda, 'Gk')
    print 'test_tcr_ss_spiking: synapses created at', datetime.now().strftime('%Y%m%d_%H%M%S')
    # Now create tables for recording Vm, [Ca2+] on soma of each cell
    # (both TCR and SS).
    for cell in chain(tcr, ss):
        vm = moose.Table('%s/vm_soma_%s' % (sim.data.path, cell.name))
        vm.stepMode = 3
        vm.connect('inputRequest', cell.soma, 'Vm')
        vm_tabs.append(vm)
        ca = moose.Table('%s/ca_soma_%s' % (sim.data.path, cell.name))
        ca.stepMode = 3
        ca.connect('inputRequest', moose.CaConc('%s/CaPool' % (cell.soma.path)), 'Ca')
        ca_tabs.append(ca)
    # Create PulseGen for stimulating each TCR cell.
    stim = [moose.PulseGen('%s/stim_%d' % (sim.model.path, ii)) for ii in range(num_tcr)]
    pulses = [1.0, 
              3.000, 3.040, 3.080, 3.120, 3.160, 3.200,
              4.0, 4.2, 4.4,
              1e9]
    width = 2e-3
    level = 2e-9
    diff = 0.0
    stim_tabs = []
    # Create stimuli that are slightly off each other by {offset}
    for st, cell in zip(stim, tcr):
        st.setCount(len(pulses)+1)
        st.level[0] = level
        st.delay[0] = pulses[0] + diff
        st.width[0] = width
        for ii, delay in enumerate(np.diff(pulses)):
            st.level[ii+1] = level
            st.delay[ii+1] = delay + diff
            st.width[ii] = width
        st.connect('outputSrc', cell.soma, 'injectMsg')
        stim_tab = moose.Table('%s/stim' % (sim.data.path))
        stim_tab.stepMode = 3
        stim_tab.connect('inputRequest', st, 'output')
        stim_tabs.append(stim_tab)
        diff += offset
        
    sim.schedule()    
    print 'test_tcr_ss_spiking: scheduling done at', datetime.now().strftime('%Y%m%d_%H%M%S')
    sim.run(.001)
    print 'test_tcr_ss_spiking: simulation done at', datetime.now().strftime('%Y%m%d_%H%M%S')
    for index, tablist in enumerate((nmda_tabs, ampa_tabs, vm_tabs, ca_tabs, stim_tabs)):
        for tab in tablist:
            ts = np.linspace(0, sim.simtime, len(tab))
            data = np.asarray(tab)
            np.savetxt('%s/%s.dat' % (datadir, tab.name), np.c_[ts, tab])            
        # pylab.legend()
        # bbox = matplotlib.transforms.Bbox.from_bounds(.1, .5, .5, .3) 
        # trans = ax.transAxes + fig.transFigure.inverted() 
        # l, b, w, h = matplotlib.transforms.TransformedBbox(bbox, trans).bounds
        # axins = fig.add_axes([l, b, w, h]) 
        # axins.plot(ts, stimdata, label='stimulus')
        # axins.set_ylim(-1e-9, 2e-9)
    # Save figures for Vm on soma of each cell
    for tab in vm_tabs:
        ts = np.linspace(0, sim.simtime, len(tab))
        plt.figure()
        plt.plot(ts, tab)
        plt.savefig('%s/%s.png' % (datadir, tab.name))
        plt.close()
    print 'test_tcr_ss_spiking: finished at', datetime.now().strftime('%Y%m%d_%H%M%S')

        

if __name__ == '__main__':
    # test_tcr_spinstell_ampa()
    test_tcr_ss_spiking(0.0)
    
# 
# test_tcr_spinstell.py ends here
