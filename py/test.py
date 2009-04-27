# test.py --- 
# 
# Filename: test.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sat Apr 18 01:08:37 2009 (+0530)
# Version: 
# Last-Updated: Mon Apr 27 19:54:20 2009 (+0530)
#           By: subhasis ray
#     Update #: 657
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
import pylab
import moose
import config
from nachans import *
from kchans import *
from cachans import CaL, CaT
from archan import AR
from capool import CaPool
from compartment import *

conductance = {'NaF': 1500.0,
               'NaF_TCR': 1500.0,
               'NaF2': 1500.0,
               'NaP': 1.5,
               'NaPF': 1.5,
               'NaPF_SS': 1.5,
               'NaPF_TCR': 1.5,
               'KDR': 1000.0,
               'KDR_FS': 1000.0,
               'CaT': 1.0,
               'CaL': 5.0,
               'KA': 300.0,
               'KA_IB': 300.0,
               'KC': 100.0,
               'KC_FAST': 100.0,
               'KM': 37.5,
               'K2': 1.0,
               'KAHP': 1.0,
               'KAHP_SLOWER': 1.0,
               'AR': 2.5}




def setup_singlecomp(channels, densities=None):
    """channels is the list of channel class names (string)"""
    if config.context.exists('test'):
        print "Model is already set up. Returning."
        container = moose.Neutral("test")
        data = moose.Neutral("data")
        return (container, data)

    container = moose.Neutral("test")
    data = moose.Neutral("data")
    tables = []
    comp = MyCompartment("comp", container)
    container.comp = comp
    comp.length = 20e-6
    comp.diameter = 15e-6
    comp.initVm = -65e-3
    comp.Em = -65e-3
    comp.setSpecificCm(9e-3)
    comp.setSpecificRm(5.0)
    comp.setSpecificRa(2.5)
    for channel in channels:
        key = None
        if type(channel) is type(''):
            key = channel
        else:
            key = channel.name
        chan = comp.insertChannel(channel, conductance[key])
        table = moose.Table(key, data)
        table.stepMode = 3
        chan.connect("Gk", table, "inputRequest")
        
    comp.insertRecorder("Vm", data)
    pulsegen = moose.PulseGen("pulsegen", container)
    pulsegen.baseLevel = 0.0
    pulsegen.firstLevel = 1e-10
    pulsegen.firstWidth = 20e-3
    pulsegen.firstDelay = 20e-3
    pulsegen.connect("outputSrc", comp, "injectMsg")
    # We cannot read comp.Im using a table because it is always set to
    # zero at the end of process method. Hence read it from the
    # pulsegen
    table = moose.Table("Inject", data)
    table.stepMode = 3
    pulsegen.connect('output', table, 'inputRequest')
    return (container, data)

class Simulation:
    """This class is a wrapper to control a whole simulation."""
    def __init__(self):
        self.model = None
        self.data = None
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


import pylab
if __name__ == "__main__":
    sim = Simulation()
    model = moose.Neutral('model')
    sim.model = model
    comp = createTestCompartment('comp', model, 
                                 20e-6, 15e-6,
                                 chan_gbar_dict={
            # 'NaF2': 750.0,
#             'NaPF_SS': 0.75,
#             'KDR_FS': 750.0,
#             'KC_FAST': 100.0,
#             'KA': 300.0,
#             'KM': 37.5,
#             'K2': 1.0,
            'KAHP_SLOWER': 1.0,
            'CaL': 5.0,
#             'CaT_A': 1.0,
#             'AR': 2.5
})
    
    comp.insertCaPool(5.2e-6/0.02e-10, 20e-3)
    comp.insertPulseGen('pulsegen', sim.model)
    data = moose.Neutral('data')
    sim.data = data
    vm_table = comp.insertRecorder('Vm', data)
    sim.schedule()
    sim.run(50e-3)
    tables = sim.dump_data('data')

##############################
    nrn_data = pylab.loadtxt('../nrn/mydata/Vm.plot')
    nrn_Vm = nrn_data[:, 1]
    nrn_t = nrn_data[:, 0]
    mus_t = pylab.array(range(len(vm_table))) * 1e-3
    pylab.plot(nrn_t, nrn_Vm, 'r-', label='nrn')
    pylab.plot(mus_t, pylab.array(vm_table)*1e3, 'g-', label='mus')
    pylab.legend()
    pylab.show()


#####################################
#     tables = test_singlecomp([], 50e-3)
#     rows = ceil(sqrt(len(tables)))
#     cols = ceil(len(tables) / float(rows))
#     for fignum in range(len(tables)):
#         pylab.subplot(rows, cols, fignum+1)
#         pylab.plot(tables[fignum], label=tables[fignum].name)
#         pylab.legend()
#     pylab.show()
# 
# test.py ends here
