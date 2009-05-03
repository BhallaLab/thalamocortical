# test.py --- 
# 
# Filename: test.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sat Apr 18 01:08:37 2009 (+0530)
# Version: 
# Last-Updated: Mon May  4 00:01:25 2009 (+0530)
#           By: subhasis ray
#     Update #: 700
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

from math import pi, ceil, floor, sqrt
import os
from datetime import datetime

import moose
import config
from nachans import *
from kchans import *
from cachans import CaL, CaT
from archan import AR
from capool import CaPool
from compartment import *
from simulation import Simulation

def createTestCompartment(name, parent, 
                          length=20e-6, diameter=15e-6, 
                          Em=-65e-3, initVm=-65e-3, CM=9e-3, RM=5.0, RA=2.5, 
                          channel_dict={}):
    """Creates a compartment with the given parameters and inserts the
    channels in channel_dict.

    The channel_dict should a dictionary containing (channel class name,
    conductance density) key-value pairs.
    """
    comp = MyCompartment(name, parent)
    comp.length = length
    comp.diameter = diameter
    comp.Em = Em
    comp.initVm = initVm
    comp.setSpecificCm(CM)
    comp.setSpecificRm(RM)
    comp.setSpecificRa(RA)
    for channel, gbar in channel_dict.items():
        channel = comp.insertChannel(channel, gbar)
        setattr(comp, channel.name, channel)
    return comp

test_conductances_ss = {
    'NaF2': 750.0,
#     'NaPF_SS': 0.75,
#     'KDR_FS': 750.0,
#     'KC_FAST': 100.0,
#     'KA': 300.0,
#     'KM': 37.5,
#     'K2': 1.0,
#     'KAHP_SLOWER': 1.0,
#     'CaL': 5.0,
#     'CaT_A': 1.0,
#     'AR': 2.5
    }


def create_testcomp(name, parent, length=20e-6, diameter=15e-6, 
                    Em=-65e-3, initVm=-65e-3, CM=9e-3, RM=5.0, RA=2.5,
                    conductance_dict={}):
    comp = MyCompartment(name, parent)
    comp.length = length
    comp.diameter = diameter
    comp.initVm = initVm
    comp.Em = Em
    comp.setSpecificCm(CM)
    comp.setSpecificRm(RM)
    comp.setSpecificRa(RA)
    for channel, density in conductance_dict.items():
        shift = None
        if channel == 'NaF2' or channel == 'NaPF_SS':
            shift = -2.5e-3
            print 'Shift set to:', shift
        print 'now :', shift
        comp.insertChannel(channel, specificGbar=density, shift=shift)
    return comp


import pylab
if __name__ == "__main__":
    sim = Simulation()
    comp = create_testcomp('comp', sim.model, 
                           conductance_dict={ 
            'NaF2': 750.0,
            'NaPF_SS': 0.75,
            'KDR_FS': 750.0,
            'KC_FAST': 100.0,
            'KA': 300.0,
            'KM': 37.5,
            'K2': 1.0,
            'KAHP_SLOWER': 1.0,
            'CaL': 5.0,
            'CaT_A': 1.0,
            'AR': 2.5
            })
    
    comp.insertPulseGen('pulsegen', sim.model)
    comp.insertCaPool(5.2e-6 / 2e-10, 20e-3) # The fortran code uses 2e-4 um depth
    ca_table = moose.Table('Ca', sim.data)
    ca_table.stepMode = 3
    comp.ca_pool.connect('Ca', ca_table, 'inputRequest')
    vm_table = comp.insertRecorder('Vm', sim.data)
    
    sim.schedule()
    sim.run(50e-3)
    tables = sim.dump_data('data')

###################################################################################
# Assumption is that a proper Vm recording is already available in
# '../nrn/mydata/Vm.plot'. This is generated by running the test.hoc
# code (properly edited) in neuron.
###################################################################################
    nrn_data = pylab.loadtxt('../nrn/mydata/Vm.plot')
    nrn_Vm = nrn_data[:, 1]
    nrn_t = nrn_data[:, 0]
    mus_t = pylab.array(range(len(vm_table)))*1e-3
    nrn_Ca = pylab.loadtxt('../nrn/mydata/Ca.plot')[:,1]
##############################
#     mus_m = pylab.array(m_table)
#     pylab.plot(mus_Ca * 1e3, mus_m)
#     pylab.plot(nrn_Ca, nrn_m)
#     pylab.show()
# ###############
    pylab.subplot(2, 1, 1, title='Vm')
    pylab.plot(nrn_t, nrn_Vm, 'rx', label='nrn')
    pylab.plot(mus_t, pylab.array(vm_table)*1e3, 'g-', label='mus')
#     pylab.legend()
    pylab.subplot(2, 1, 2, title='[Ca2+]')
    pylab.plot(nrn_t, nrn_Ca, 'rx', label='nrn')
    pylab.plot(mus_t, pylab.array(ca_table) * 1e3, 'g-', label='mus')
    pylab.legend()

    pylab.show()


# 
# test.py ends here
