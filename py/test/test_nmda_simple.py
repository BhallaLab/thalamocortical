# test_nmda.py --- 
# 
# Filename: test_nmda.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Mar 22 16:58:57 2010 (+0530)
# Version: 
# Last-Updated: Sat Dec 29 18:27:12 2012 (+0530)
#           By: subha
#     Update #: 197
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This code is for testing the NMDA synapse in MOOSE - the comparison is
# nmda.mod in neuron version of Traub et al 2005 model. The test file
# for neuron is test_nmda.hoc
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
import pylab
import numpy as np
import moose

def testNMDAChan(simtime=1.0, simdt=1e-5, plotdt=1e-5):
    context = moose.PyMooseBase.getContext()
    container = moose.Neutral('test_NMDA')

    soma_b = moose.Compartment('B', container)
    soma_b.Rm = 5.3e9 # GM = 2e-5 S/cm^2
    soma_b.Cm = 8.4823001646924426e-12 # CM = 0.9 uF/cm^2
    soma_b.Em = -65e-3 
    soma_b.initVm = -65e-3
    soma_b.Ra = 282942.12 # RA = 250 Ohm-cm
    
    nmda = moose.NMDAChan('nmda', container)
    nmda.tau2 = 5e-3
    nmda.tau1 = 130e-3
    nmda.Gbar = 1e-9
    nmda.saturation = 1e10
    nmda.connect('channel', soma_b, 'channel')
    spikegen = moose.SpikeGen('spike', container)
    spikegen.threshold = 0.5
    spikegen.connect('event', nmda, 'synapse')
    spikegen.refractT = 0.0
    nmda.delay[0] = 1e-3
    nmda.weight[0] = 1.0
    
    pulsegen = moose.PulseGen('pulse', container)
    pulsegen.setCount(3)
    pulsegen.level[0] = 1.0
    pulsegen.delay[0] = 10e-3
    pulsegen.width[0] = 1e-3
    pulsegen.level[1] = 1.0
    pulsegen.delay[1] = 2e-3
    pulsegen.width[1] = 1e-3    
    pulsegen.delay[2] = 1e9
    pulsegen.connect('outputSrc', spikegen, 'Vm')

    data = moose.Neutral('data')
    vmB = moose.Table('Vm_B', data)
    vmB.stepMode = 3
    vmB.connect('inputRequest', soma_b, 'Vm')
    pulse = moose.Table('pulse', data)
    pulse.stepMode = 3
    pulse.connect('inputRequest', pulsegen, 'output')
    gNMDA = moose.Table('G', data)
    gNMDA.stepMode = 3
    gNMDA.connect('inputRequest', nmda, 'Gk')

    context.setClock(0, simdt)
    context.setClock(1, simdt)
    context.setClock(2, simdt)
    context.setClock(3, plotdt)
    context.reset()
    context.step(simtime)
    # gNMDA.dumpFile('gNMDA.dat', False)
    # vmA.dumpFile('Va.dat', False)
    # vmB.dumpFile('Vb.dat', False)
    ts = np.linspace(0, simtime, len(gNMDA))
    pylab.plot(ts, pulse)
    pylab.plot(ts, np.asarray(gNMDA) * 1e9, label='gNMDA')
    pylab.show()
    np.savetxt('../data/two_comp_nmda.plot', np.transpose(np.vstack((ts, vmB, gNMDA))))

if __name__ == '__main__':
    testNMDAChan()

    
    

# 
# test_nmda.py ends here
