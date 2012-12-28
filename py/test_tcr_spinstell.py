# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Fri Dec 28 10:52:22 2012 (+0530)
#           By: subha
#     Update #: 173
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

def test_tcr_ss_spiking():
    netdata = TraubFullNetData()
    sim = Simulation('tcr_ss')
    tcr_idx = netdata.celltype.index('TCR')
    ss_idx = netdata.celltype.index('SpinyStellate')
    num_ss, num_tcr = 1, 20
    ss = [SpinyStellate(SpinyStellate.prototype, '%s/SS_%d' % (sim.model.path, idx)) for idx in range(num_ss)]
    tcr = [TCR(TCR.prototype, '%s/TCR_%d' % (sim.model.path, idx)) for idx in range(num_tcr)]
    pre_per_post = netdata.pre_post_ratio[tcr_idx][ss_idx]
    for cell in ss:
        post_comp_list = [cell.comp[ii] for ii in random.sample(netdata.allowed_comps[tcr_idx][ss_idx], pre_per_post)]
        nmda_tabs = []
        ampa_tabs = []
        vm_tabs = []
        ca_tabs = []
        for precell, postcomp in izip(tcr, cycle(post_comp_list)):
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
    for cell in chain(tcr, ss):
        vm = moose.Table('%s/vm_soma_%s' % (sim.data.path, cell.name))
        vm.stepMode = 3
        vm.connect('inputRequest', cell.soma, 'Vm')
        vm_tabs.append(vm)
        ca = moose.Table('%s/ca_soma_%s' % (sim.data.path, cell.name))
        ca.stepMode = 3
        ca.connect('inputRequest', moose.CaConc('%s/CaPool' % (cell.soma.path)), 'Ca')
        ca_tabs.append(ca)
    stim = moose.PulseGen('%s/stim')
    pulsecount = 6
    delay = 25e-3
    width = 2e-3
    level = 0.3e-9
    stim.setCount(7)
    for ii in range(pulsecount):
        stim.level[ii] = level
        stim.delay[ii] = delay
        stim.width[ii] = width
    stim.delay[0] = 1.0
    sim.schedule()
    sim.run(5.0)        
    for tab in chain(nmda_tabs, ampa_tabs, vm_tabs, ca_tabs):
        ts = np.linspace(0, sim.simtime, len(tab))
        np.savetxt('%s.dat' % (tab.name), np.c_[ts, tab])
    print 'Finished'
        

def test_tcr_spinstell_ampa():
    netdata = TraubFullNetData()
    sim = Simulation('tcr_spinstell_synapse')
    tcr_index = netdata.celltype.index('TCR')
    spinstell_index = netdata.celltype.index('SpinyStellate')
    tcr = TCR(TCR.prototype, sim.model.path + '/TCR')
    spinstell = SpinyStellate(SpinyStellate.prototype, sim.model.path + '/SpinyStellate')
    precomp = tcr.comp[TCR.presyn]
    postcomp = spinstell.comp[31] # 5 is among the allowed post synaptic compartments in spiny stellate cell
    tau_ampa = netdata.tau_ampa[tcr_index][spinstell_index]
    synchan = precomp.makeSynapse(postcomp, 
                                  name='ampa_from_TCR', 
                                  classname='SynChan', 
                                  Ek=0.0,
                                  Gbar=netdata.g_ampa_baseline[tcr_index][spinstell_index] * tau_ampa*1e3/pylab.e,
                                  tau1=tau_ampa,
                                  tau2=tau_ampa,
                                  delay = synapse.SYNAPTIC_DELAY_THALAMOCORTICAL
                                  )
    stim = tcr.soma.insertPulseGen('stimulus', sim.model, firstLevel=1e-9, firstDelay=200e-3, firstWidth=2e-3)
    tcr_soma_tab = tcr.soma.insertRecorder('stim', 'Vm', sim.data)
    ss_soma_tab = spinstell.soma.insertRecorder('ss_soma', 'Vm', sim.data)
    ss_dend_tab = postcomp.insertRecorder('ss_dend', 'Vm', sim.data)
    gk_ampa_tab = moose.Table('gk_ss', sim.data)
    gk_ampa_tab.stepMode = 3
    print 'Connected Gk', gk_ampa_tab.connect('inputRequest', synchan, 'Gk')
    sim.schedule()
    sim.run(1.0)
    pylab.plot(np.linspace(0, 1.0, len(tcr_soma_tab)), tcr_soma_tab, label='tcr_soma')
    pylab.plot(np.linspace(0, 1.0, len(tcr_soma_tab)), ss_soma_tab, label='ss_soma')
    pylab.plot(np.linspace(0, 1.0, len(tcr_soma_tab)), ss_dend_tab, label='ss_dend')
    pylab.plot(np.linspace(0, 1.0, len(gk_ampa_tab)), np.array(gk_ampa_tab) * 1e9, label='gk_ampa_spinstell (nS)')
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    # test_tcr_spinstell_ampa()
    test_tcr_ss_spiking()
    
# 
# test_tcr_spinstell.py ends here
