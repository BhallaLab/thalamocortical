# test_tcr_hsolve.py --- 
# 
# Filename: test_tcr_hsolve.py
# Description: 
# Author: 
# Maintainer: 
# Created: Thu Jul 11 11:45:52 2013 (+0530)
# Version: 
# Last-Updated: Thu Jul 11 15:23:11 2013 (+0530)
#           By: subha
#     Update #: 126
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

"""This script tests the TCR cell with a current injection: for comparison with new moose hsolve."""

import os
from unittest import TestCase
import unittest
from datetime import datetime
import numpy as np
import moose
import config
from tcr import TCR

simtime = 100e-3
simdt = 1e-8
plotdt = 0.25e-3

pulsearray = [[50e-3, 20e-3, 1e-9],
              [1e9, 0, 0]]


class TestTCRHSolveEE(TestCase):
    def setUp(self):
        self.model = moose.Neutral('/model')
        self.data = moose.Neutral('/data')
        self.tcr = TCR(TCR.prototype, self.model.path + '/TCR')
        self.stimulus = self.tcr.soma.insertPulseGen('pulse', self.model)        
        self.soma_vm = self.tcr.soma.insertRecorder('somaVm', 'Vm', self.data)
        self.presyn_vm = self.tcr.comp[self.tcr.presyn].insertRecorder('presynVm', 'Vm', self.data)
        self.stimtab = moose.Table('%s/injectionCurrent' % (self.data.path))
        self.stimtab.stepMode = 3
        self.stimtab.connect('inputRequest', self.stimulus, 'output')

    def setstimulus(self, pulsearray):
        self.stimulus.count = len(pulsearray)
        for ii, (delay, width, level) in enumerate(pulsearray):
            self.stimulus.delay[ii] = delay
            self.stimulus.width[ii] = width
            self.stimulus.level[ii] = level
        
    def schedule(self, simdt, plotdt, solver):
        config.LOGGER.info('Scheduling: simdt=%g, plotdt=%g, solver=%s' % (simdt, plotdt, solver))
        self.simdt = simdt
        self.plotdt = plotdt
        self.solver = solver
        self.tcr.method = solver
        config.context.setClock(0, self.simdt)
        config.context.setClock(1, self.simdt)
        config.context.setClock(2, self.simdt)
        config.context.setClock(3, self.simdt)
        config.context.setClock(4, self.plotdt)
        config.context.useClock(0, '%s/##[TYPE=Compartment]' % (self.model.path), 'init')
        config.context.useClock(1, '%s/##[TYPE=Compartment]' % (self.model.path), 'process')
        config.context.useClock(3, '%s/##[TYPE!=Compartment]' % (self.model.path), 'process')
        config.context.useClock(4, '%s/##[TYPE=Table]' % (self.data.path))
        if self.solver == 'hsolve':
            config.context.useClock(2, '%s/##[TYPE=HSolve]' % (self.model.path))
    
    def runsim(self, simtime, stepsize):
        self.simtime = simtime
        config.context.reset()
        config.LOGGER.info('Running simulation: simtime=%g, stepsize=%g' % (simtime, stepsize))
        t = 0
        tleft = simtime
        while tleft > stepsize:
            ts = datetime.now()
            config.context.step(stepsize)
            te = datetime.now()
            td = te - ts
            t += stepsize
            tleft -= stepsize
            config.LOGGER.info('Run till t=%g, left: %g, time to run %g s simulationwith simdt=%g: %g' % (t, tleft, stepsize, self.simdt, td.days * 86400 + td.seconds + td.microseconds*1e-6))
        if tleft > 0:
            ts = datetime.now()
            config.context.step(tleft)
            te = datetime.now()
            td = te - ts
            config.LOGGER.info('Run till t=%g, time to run %g s simulation with simdt=%g and solver=%s: %g' % (self.simtime, tleft, self.simdt, self.solver, td.days * 86400 + td.seconds + td.microseconds*1e-6))

    def savedata(self):
        for tid in self.data.children():
            tab = moose.Table(tid)
            ts = np.linspace(0, self.simtime, len(tab))
            data = np.vstack((ts, np.asarray(tab))).transpose()
            fname = '%s_%s_%s_%s.dat' % (self.tcr.name, tab.name, self.solver, config.filename_suffix) 
            path = os.path.join(config.data_dir, fname)
            np.savetxt(path, data)
            config.LOGGER.info('Saved table %s in file %s' % (tab.name, path))

    def testHSolve(self):
        self.setstimulus(pulsearray)
        self.schedule(simdt, plotdt, 'hsolve')
        self.runsim(simtime, simtime/10)
        self.savedata()

    def testEE(self):
        self.setstimulus(pulsearray)
        self.schedule(simdt, plotdt, 'ee')
        self.runsim(simtime, simtime/10)
        self.savedata()
        

if __name__ == '__main__':
    unittest.main()

# 
# test_tcr_hsolve.py ends here
