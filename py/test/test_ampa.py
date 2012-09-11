# test_ampa.py --- 
# 
# Filename: test_ampa.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Mar 22 16:58:57 2010 (+0530)
# Version: 
# Last-Updated: Tue Sep 11 17:42:19 2012 (+0530)
#           By: subha
#     Update #: 191
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This code is for testing the AMPA synapse in MOOSE - the comparison is
# ampa.mod in neuron version of Traub et al 2005 model. The test file
# for neuron is test_ampa.hoc
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

import math
import numpy as np
import pylab
import moose

def testAMPAChan(simtime=1000e-3, simdt=1e-5, plotdt=1e-5):
    context = moose.PyMooseBase.getContext()
    container = moose.Neutral('test_AMPA')
    soma_a = moose.Compartment('A', container)
    soma_a.Rm = 5.3e9 # GM = 2e-5 S/cm^2
    soma_a.Cm = 8.4823001646924426e-12 # CM = 0.9 uF/cm^2
    soma_a.Em = -65e-3
    soma_a.initVm = -65e-3
    soma_a.Ra = 282942.12
    
    soma_b = moose.Compartment('B', container)
    soma_b.Rm = 5.3e9 # GM = 2e-5 S/cm^2
    soma_b.Cm = 8.4823001646924426e-12 # CM = 0.9 uF/cm^2
    soma_b.Em = -65e-3 
    soma_b.initVm = -65e-3
    soma_b.Ra = 282942.12 # RA = 250 Ohm-cm
    
    ampa = moose.SynChan('ampa', container)
    ampa.tau2 = 2e-3 
    ampa.tau1 = ampa.tau2
    ampa.connect('channel', soma_b, 'channel')
    
    spikegen = moose.SpikeGen('spike', container)
    spikegen.threshold = 0.0
    spikegen.edgeTriggered = 1
    # spikegen.absRefractT = 1
    
    soma_a.connect('VmSrc', spikegen, 'Vm')
    spikegen.connect('event', ampa, 'synapse')

    ampa.delay[0] = 5e-3
    # Traub et al use the equation gAMPA = c * t * exp(-t/tau), making
    # the peak conductance c * tau / e. This means c will have unit of
    # conductance/time. Hence scaling from 0.5e-3 uS/ms to SI.
    ampa.Gbar = 0.5e-3 * 1e-6 * ampa.tau1*1e3 / math.e 
    ampa.weight[0] = 1.0

    pulsegen = moose.PulseGen('pulse', container)
    pulsegen.firstLevel = 0.1e-9
    pulsegen.firstDelay = 20e-3
    pulsegen.firstWidth = 20e-3
    pulsegen.secondLevel = 0.0
    pulsegen.secondDelay = 1e9
    pulsegen.connect('outputSrc', soma_a, 'injectMsg')

    data = moose.Neutral('data')
    vmB = moose.Table('Vm_B', data)
    vmB.stepMode = 3
    vmB.connect('inputRequest', soma_b, 'Vm')

    vmA = moose.Table('Vm_A', data)
    vmA.stepMode = 3
    vmA.connect('inputRequest', soma_a, 'Vm')

    gAMPA = moose.Table('G_ampa', data)
    gAMPA.stepMode = 3
    gAMPA.connect('inputRequest', ampa, 'Gk')

    context.setClock(0, simdt)
    context.setClock(1, simdt)
    context.reset()
    context.step(simtime)
    ts = np.linspace(0, simtime, len(gAMPA))
    np.savetxt('../data/two_comp_ampa.plot', np.transpose(np.vstack((ts, vmA, vmB, gAMPA))))
    # gAMPA.dumpFile('gAMPA.dat', False)
    # vmA.dumpFile('Va.dat', False)
    # vmB.dumpFile('Vb.dat', False)
    pylab.plot(ts, gAMPA)
    pylab.show()

if __name__ == '__main__':
    testAMPAChan()

    
    

# 
# test_ampa.py ends here
