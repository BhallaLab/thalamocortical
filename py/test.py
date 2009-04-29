# test.py --- 
# 
# Filename: test.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sat Apr 18 01:08:37 2009 (+0530)
# Version: 
# Last-Updated: Thu Apr 30 02:25:30 2009 (+0530)
#           By: subhasis ray
#     Update #: 653
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

class Simulation:
    """This class is a wrapper to control a whole simulation."""
    def __init__(self):
        self.model = moose.Neutral('model')
        self.data = moose.Neutral('data')
        self.start_t = None
        self.end_t = None

    def schedule(self):
        config.context.setClock(0, config.simdt)
        config.context.setClock(1, config.simdt)
        config.context.setClock(2, config.plotdt)

    def run(self, time):
        config.context.reset()
        self.start_t = datetime.now()
        config.context.step(float(time))
        self.end_t = datetime.now()

    def dump_data(self, directory, time_stamp=False):
        """Save the data in directory. It creates a subdirectory with
        the date of start_t. The files are prefixed with start time in
        HH.MM.SS_ format if time_stamp is True."""
        path = directory
        tables = []
        if not os.access(directory, os.W_OK):
            print 'data directory:', directory, 'is not writable'
            return
        else:
            path = directory + '/' + self.start_t.strftime('%Y_%m_%d') + '/'
            if not os.access(path, os.F_OK):
                os.mkdir(path)
        if time_stamp:
            path = path + self.start_t.strftime('%H.%M.%S') + '_'
        for table_id in self.data.children():
            table = moose.Table(table_id)
            tables.append(table)
            file_path = path + table.name + '.plot'
            table.dumpFile(file_path)
            print 'Dumped data in ', file_path
        return tables

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
#     'NaF2': 750.0,
#     'NaPF_SS': 0.75,
#     'KDR_FS': 750.0,
#     'KC_FAST': 100.0,
#     'KA': 300.0,
#     'KM': 37.5,
#     'K2': 1.0,
#     'KAHP_SLOWER': 1.0,
    'CaL': 5.0,
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
        comp.insertChannel(channel, density)
    return comp


import pylab
if __name__ == "__main__":
    sim = Simulation()
    comp = create_testcomp('comp', sim.model, 
                           conductance_dict={ 'NaF2': 750.0,
                                              'NaPF_SS': 0.75,
                                              'KDR_FS': 750.0,
                                              'KC_FAST': 100.0,
                                              'KA': 300.0,
                                              'KM': 37.5,
                                              'K2': 1.0,
                                              'KAHP_SLOWER': 1.0,
                                              'CaL': 5.0,
                                              'CaT_A': 1.0,
                                              'AR': 2.5})
    
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
    pylab.legend()
    pylab.subplot(2, 1, 2, title='[Ca2+]')
    pylab.plot(nrn_t, nrn_Ca, 'rx', label='nrn')
    pylab.plot(mus_t, pylab.array(ca_table) * 1e3, 'g-', label='mus')
    pylab.legend()
#     pylab.legend()
    pylab.show()


# 
# test.py ends here
