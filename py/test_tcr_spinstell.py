# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Wed Jan 18 21:56:10 2012 (+0530)
#           By: subha
#     Update #: 73
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

import numpy
import pylab

import moose

from simulation import Simulation
from trbnetdata import TraubFullNetData
from spinystellate import SpinyStellate
from tcr import TCR
import synapse

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
    pylab.plot(numpy.linspace(0, 1.0, len(tcr_soma_tab)), tcr_soma_tab, label='tcr_soma')
    pylab.plot(numpy.linspace(0, 1.0, len(tcr_soma_tab)), ss_soma_tab, label='ss_soma')
    pylab.plot(numpy.linspace(0, 1.0, len(tcr_soma_tab)), ss_dend_tab, label='ss_dend')
    pylab.plot(numpy.linspace(0, 1.0, len(gk_ampa_tab)), numpy.array(gk_ampa_tab) * 1e9, label='gk_ampa_spinstell (nS)')
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    test_tcr_spinstell_ampa()
    
    
# 
# test_tcr_spinstell.py ends here
