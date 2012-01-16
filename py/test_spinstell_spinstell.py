# test_tcr_spinstell.py --- 
# 
# Filename: test_tcr_spinstell.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jan 16 09:50:05 2012 (+0530)
# Version: 
# Last-Updated: Mon Jan 16 11:31:57 2012 (+0530)
#           By: subha
#     Update #: 74
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
import synapse

def test_spinstell_spinstell_ampa():
    netdata = TraubFullNetData()
    sim = Simulation('spinstell_spinstell_synapse')
    spinstell_index = netdata.celltype.index('SpinyStellate')
    pre = SpinyStellate(SpinyStellate.prototype, sim.model.path + '/SpinyStellate1')
    spinstell = SpinyStellate(SpinyStellate.prototype, sim.model.path + '/SpinyStellate2')
    precomp = pre.comp[SpinyStellate.presyn]
    postcomp = spinstell.comp[5] # 5 is among the allowed post synaptic compartments in spiny stellate cell
    synchan = precomp.makeSynapse(postcomp, 
                                  name='ampa_from_SPINSTELL', 
                                  classname='SynChan', 
                                  Ek=0.0,
                                  Gbar=netdata.g_ampa_baseline[spinstell_index][spinstell_index],
                                  tau1=netdata.tau_ampa[spinstell_index][spinstell_index],
                                  tau2=netdata.tau_ampa[spinstell_index][spinstell_index],
                                  delay = synapse.SYNAPTIC_DELAY_DEFAULT
                                  )
    stim = pre.soma.insertPulseGen('stimulus', sim.model, firstLevel=1e-9, firstDelay=200e-3, firstWidth=2e-3)
    pre_soma_tab = pre.soma.insertRecorder('stim', 'Vm', sim.data)
    ss_soma_tab = spinstell.soma.insertRecorder('post_soma', 'Vm', sim.data)
    ss_dend_tab = postcomp.insertRecorder('post_dend', 'Vm', sim.data)
    sim.schedule()
    sim.run(1.0)
    pylab.plot(numpy.linspace(0, 1.0, len(pre_soma_tab)), pre_soma_tab, label='pre_soma')
    pylab.plot(numpy.linspace(0, 1.0, len(ss_soma_tab)), ss_soma_tab, label='ss_soma')
    pylab.plot(numpy.linspace(0, 1.0, len(ss_dend_tab)), ss_dend_tab, label='ss_dend')
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    test_spinstell_spinstell_ampa()
    
    
# 
# test_tcr_spinstell.py ends here
