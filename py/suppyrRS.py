# suppyrrs.py --- 
# 
# Filename: suppyrrs.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Aug  7 13:59:30 2009 (+0530)
# Version: 
# Last-Updated: Mon Apr 26 17:55:49 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 661
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: Superficial Regular Spiking Pyramidal Cells of layer 2/3
# From Traub et all, 2005
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

from datetime import datetime
import config
import trbutil
import moose
from cell import *
from capool import CaPool

class SupPyrRS(TraubCell):
    """Superficial Pyramidal Regula Spiking cell."""
    chan_params = {
        'ENa': 50e-3,
        'EK': -95e-3,
        'ECa': 125e-3,
        'EAR': -35e-3,
        'EGABA': -81e-3,
        'TauCa': 20e-3
        }
    ca_dep_chans = ['KAHP', 'KC']
    num_comp = 74
    presyn = 72
    proto_file = "SupPyrRS.p"
    prototype = TraubCell.read_proto(proto_file, "SupPyrRS", chan_params)
    def __init__(self, *args):
        start = datetime.now()
	TraubCell.__init__(self, *args)
        soma_ca_pool = moose.CaConc(self.soma.path + '/CaPool')
        soma_ca_pool.tau = 100e-3
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('created cell in: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        # self._setup_channels()
	
    def _topology(self):
        raise Exception, 'Deprecated'
        self.presyn = 72
	self.level[1].add(self.comp[1])
	for ii in range(2,14):
	    self.level[2].add(self.comp[ii])
	for ii in range(14, 26):
	    self.level[3].add(self.comp[ii])
	for ii in range(26, 38):
	    self.level[4].add(self.comp[ii])
	self.level[5].add(self.comp[38])
	self.level[6].add(self.comp[39])
	self.level[7].add(self.comp[40])
	self.level[8].add(self.comp[41])
	self.level[8].add(self.comp[42])
	self.level[9].add(self.comp[43])
	self.level[9].add(self.comp[44])
	for ii in range(45, 53):
	    self.level[10].add(self.comp[ii])
	for ii in range(53, 61):
	    self.level[11].add(self.comp[ii])
	for ii in range(61, 69):
	    self.level[12].add(self.comp[ii])
	for ii in range(69, 75):
	    self.level[0].add(self.comp[ii])
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'
        for comp in self.comp[1:]:
            comp.Em = -70e-3

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'
        for ii in range(SupPyrRS.num_comp):
            comp = self.comp[ii+1]
            ca_pool = None
            ca_dep_chans = []
            ca_chans = []
            for child in comp.children():
                obj = moose.Neutral(child)
                if obj.name == 'CaPool':
                    ca_pool = moose.CaConc(child)
                    ca_pool.B = ca_pool.B * 1e3
                    ca_pool.tau = 1e-3/0.05
                else:
                    obj_class = obj.className
                    if obj_class == "HHChannel":
                        obj = moose.HHChannel(child)
                        pyclass = eval(obj.name)                            

                        if issubclass(pyclass, KChannel):
                            obj.Ek = -95e-3
                            if issubclass(pyclass, KCaChannel):
                                ca_dep_chans.append(obj)
                        elif issubclass(pyclass, NaChannel):
                            obj.Ek = 50e-3
                        elif issubclass(pyclass, CaChannel):
                            obj.Ek = 125e-3
                            if issubclass(pyclass, CaL):
                                ca_chans.append(obj)
                        elif issubclass(pyclass, AR):
                            obj.Ek = -35e-3
                # if ca_pool: # Setup connections for CaPool : from CaL, to KAHP and KC
                #     for channel in ca_chans:
                #         channel.connect('IkSrc', ca_pool, 'current')

                #     for channel in ca_dep_chans:
                #         channel.useConcentration = 1
                #         ca_pool.connect("concSrc", channel, "concen")



        obj = moose.CaConc(self.soma.path + '/CaPool')
        obj.tau = 1e-3 / 0.01
        config.LOGGER.debug(obj.path + 'set tau to %g' % (obj.tau))

    @classmethod
    def test_single_cell(cls):
        """Simulates a single superficial pyramidal regular spiking
        cell and plots the Vm and [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = SupPyrRS(SupPyrRS.prototype, sim.model.path + "/SupPyrRS")
        config.LOGGER.info('Created cell: %s' % (mycell.path))
        vm_table = mycell.comp[SupPyrRS.presyn].insertRecorder('Vm_suppyrrs', 'Vm', sim.data)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=50e-3, firstWidth=50e-3)
#         pulsegen1 = mycell.soma.insertPulseGen('pulsegen1', sim.model, firstLevel=3e-7, firstDelay=150e-3, firstWidth=10e-3)

        sim.schedule()
        if mycell.has_cycle():
            config.LOGGER.warning("WARNING!! CYCLE PRESENT IN CICRUIT.")
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        config.BENCHMARK_LOGGER.info('simulation time: %g' % (delta.seconds + 1e-6 * delta.microseconds))
        # sim.dump_data('data')
        # mycell.dump_cell('suppyrrs.txt')
        
        mus_vm = pylab.array(vm_table) * 1e3
        nrn_vm = pylab.loadtxt('../nrn/mydata/Vm_suppyrrs.plot')
        nrn_t = nrn_vm[:, 0]
        mus_t = linspace(0, nrn_t[-1], len(mus_vm))
        nrn_vm = nrn_vm[:, 1]
        nrn_ca = pylab.loadtxt('../nrn/mydata/Ca_suppyrrs.plot')
        nrn_ca = nrn_ca[:,1]
        pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
        pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
        # if ca_table:
        #     ca_array = pylab.array(ca_table)
        #     pylab.plot(nrn_t, -nrn_ca, 'r-', label='nrn (-)ca')
        #     pylab.plot(mus_t, -ca_array, 'b-.', label='mus (-)ca')
        #     print pylab.amax(ca_table)
        pylab.legend()
        pylab.show()
        

import unittest
class SupPyrRSTestCase(unittest.TestCase):
    def setUp(self):
        self.sim = Simulation('SupPyrRS')
        path = self.sim.model.path + '/TestSupPyrRS'
        config.LOGGER.debug('Creating cell %s' % path)
        TraubCell.adjust_chanlib(SupPyrRS.chan_params)
        self.cell = SupPyrRS(path, SupPyrRS.proto_file)
        config.LOGGER.debug('Cell created')
        self.sim.schedule()

    
    def test_compartment_count(self):
        for comp_no in range(SupPyrRS.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))

    def test_initVm(self):
        for comp_no in range(SupPyrRS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -65e-3)

    def test_Em(self):
        for comp_no in range(SupPyrRS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].Em, -70e-3)

    def test_Ca_connections(self):
        for comp_no in range(SupPyrRS.num_comp):
            ca_path = self.cell.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                config.LOGGER.debug('%s : No CaPool' % (self.cell.comp[comp_no + 1].path))
                continue
            caPool = moose.CaConc(ca_path)
            for chan in SupPyrRS.ca_dep_chans:
                chan_path = self.cell.comp[comp_no + 1].path + '/' + chan
                if not config.context.exists(chan_path):
                    continue
                chan_obj = moose.HHChannel(chan_path)
                self.assertTrue(len(chan_obj.neighbours('concen')) > 0)
                
            sources = caPool.neighbours('current')
            self.failIfEqual(len(sources), 0)
            for chan in sources:
                self.assertTrue(chan.path().endswith('CaL'))
                    
    def test_reversal_potentials(self):
        for num in range(SupPyrRS.num_comp):
            comp = self.cell.comp[num + 1]
            for chan_id in comp.neighbours('channel'):
                chan = moose.HHChannel(chan_id)
                chan_class = eval(chan.name)
                key = None
                if issubclass(chan_class, NaChannel):
                    key = 'ENa'
                elif issubclass(chan_class, KChannel):
                    key = 'EK'
                elif issubclass(chan_class, CaChannel):
                    key = 'ECa'
                elif issubclass(chan_class, AR):
                    key = 'EAR'
                else:
                    pass
                self.assertAlmostEqual(chan.Ek, SupPyrRS.chan_params[key])
    
# test main --
from simulation import Simulation
import pylab
from subprocess import call
if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_suppyrRS.hoc'], cwd='../nrn')
    # unittest.main()
    SupPyrRS.test_single_cell()
    

# 
# suppyrrs.py ends here
