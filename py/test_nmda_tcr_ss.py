# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Sat Dec 29 15:01:43 2012 (+0530)
#           By: subha
#     Update #: 204
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


def test_tcr_spinstell_nmda():
    netdata = TraubFullNetData()
    sim = Simulation('tcr_spinstell_synapse')
    tcr_index = netdata.celltype.index('TCR')
    spinstell_index = netdata.celltype.index('SpinyStellate')
    tcr = TCR(TCR.prototype, sim.model.path + '/TCR')
    spinstell = SpinyStellate(SpinyStellate.prototype, sim.model.path + '/SpinyStellate')
    precomp = tcr.comp[TCR.presyn]
    postcomp = spinstell.comp[31] # 5 is among the allowed post synaptic compartments in spiny stellate cell
    tau_nmda = netdata.tau_nmda[tcr_index][spinstell_index]
    synchan = precomp.makeSynapse(postcomp, 
                                  name='nmda_from_TCR', 
                                  classname='NMDAChan', 
                                  Ek=0.0,
                                  Gbar=netdata.g_nmda_baseline[tcr_index][spinstell_index] * tau_nmda*1e3/pylab.e,
                                  tau1=tau_nmda,
                                  tau2=5e-3,
                                  delay = synapse.SYNAPTIC_DELAY_THALAMOCORTICAL
                                  )
    synchan.MgConc = 1.5
    stim = tcr.soma.insertPulseGen('stimulus', sim.model, firstLevel=1e-9, firstDelay=200-3, firstWidth=2e-3)
    tcr_soma_tab = tcr.soma.insertRecorder('stim', 'Vm', sim.data)
    ss_soma_tab = spinstell.soma.insertRecorder('ss_soma', 'Vm', sim.data)
    ss_dend_tab = postcomp.insertRecorder('ss_dend', 'Vm', sim.data)
    gk_nmda_tab = moose.Table('gk_ss', sim.data)
    gk_nmda_tab.stepMode = 3
    print 'Connected Gk', gk_nmda_tab.connect('inputRequest', synchan, 'Gk')
    simtime = 5.0
    sim.schedule()
    sim.run(simtime)
    pylab.plot(np.linspace(0, simtime, len(tcr_soma_tab)), tcr_soma_tab, label='tcr_soma')
    pylab.plot(np.linspace(0, simtime, len(tcr_soma_tab)), ss_soma_tab, label='ss_soma')
    pylab.plot(np.linspace(0, simtime, len(tcr_soma_tab)), ss_dend_tab, label='ss_dend')
    pylab.plot(np.linspace(0, simtime, len(gk_nmda_tab)), np.array(gk_nmda_tab) * 1e9, label='gk_nmda_spinstell (nS)')
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    # test_tcr_spinstell_ampa()
    test_tcr_spinstell_nmda()
    for ch in moose.context.getWildcardList('/##[ISA=SynChan]', True):
        print ch.path()
    # test_tcr_ss_spiking()
    
# 
# test_tcr_spinstell.py ends here
