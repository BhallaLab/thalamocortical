# test.py --- 
# 
# Filename: test.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sat Apr 18 01:08:37 2009 (+0530)
# Version: 
# Last-Updated: Tue Apr 21 19:17:04 2009 (+0530)
#           By: subhasis ray
#     Update #: 344
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
from cachans import *

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
               'KC': 100.0,
               'KM': 37.5,
               'K2': 1.0,
               'KAHP': 1.0,
               'AR': 2.5}


class MyCompartment(moose.Compartment):
    def __init__(self, *args):
        moose.Compartment.__init__(self, *args)
        self.channels = []
        self._xarea = None
        self._sarea = None

    def setSpecificRm(self, RM):
        self.Rm = RM / self.sarea()

    def setSpecificRa(self, RA):
        self.Ra = RA * self.length / self.xarea()

    def setSpecificCm(self, CM):
        self.Cm = CM * self.sarea()

    def xarea(self):
        if self._xarea is None:
            self._xarea = pi * self.diameter * self.diameter
        return self._xarea

    def sarea(self):
        if self._sarea is None:
            self._sarea = pi * self.length * self.diameter
        return self._sarea

    def insertChannel(self, channel, specificGbar=None, Ek=None):
        """Insert a channel setting its gbar as membrane_area *
        specificGbar and reversal potential to Ek.
        
        This method expects either a valid channel class name or an
        existing channel object. If specificGbar is given, the Gbar is
        set to specificGbar * surface-area of the compartment. If Ek
        is given, the channel's Ek is set to this value.
        """
        if type(channel) is type(''): # if it is a class name, create the channel as a child with the same name as the class name
            
            chan_class = eval(channel)
            chan = chan_class(channel, self)
        elif type(channel) is moose.HHChannel:
            chan = channel
        else:
            print "ERROR: unknown object passed as channel: ", channel
        if specificGbar is not None:
            chan.Gbar = specificGbar * self.sarea()
        if Ek is not None:
            chan.Ek = Ek
        self.channels.append(chan)
        self.connect("channel", chan, "channel")
        return chan

    def insertRecorder(self, field_name, data_container):
        """Creates a table for recording a field under data_container"""
        table = moose.Table(field_name, data_container)
        table.stepMode = 3
        self.connect(field_name, table, "inputRequest")
        return table


def setup_singlecomp(channels):
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
    sim.model, sim.data, = setup_singlecomp(['KM'])
    sim.schedule()
    sim.run(50e-3)
    tables = sim.dump_data('data')
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
