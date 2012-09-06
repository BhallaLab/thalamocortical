# suppyrrs.py --- 
# 
# Filename: suppyrrs.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Aug  7 13:59:30 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 11:07:44 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 791
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
    # level maps level number to the set of compartments belonging to it
    level = TraubCell.readlevels("SupPyrRS.levels")
    # depth stores a map between level number and the depth of the compartments.
    depth = {
        1: 850.0 * 1e-6,
        2: 885.0 * 1e-6,
        3: 920.0 * 1e-6,
        4: 955.0 * 1e-6,
        5: 825.0 * 1e-6,
        6: 775.0 * 1e-6,
        7: 725.0 * 1e-6,
        8: 690.0 * 1e-6,
        9: 655.0 * 1e-6,
        10: 620.0 * 1e-6,
        11: 585.0 * 1e-6,
        12: 550.0 * 1e-6
        }
    prototype = TraubCell.read_proto(proto_file, "SupPyrRS", level_dict=level, depth_dict=depth, params=chan_params)
    
    def __init__(self, *args):
        # start = datetime.now()
        TraubCell.__init__(self, *args)
        soma_ca_pool = moose.CaConc(self.soma.path + '/CaPool')
        soma_ca_pool.tau = 100e-3
        # end = datetime.now()
        # delta = end - start
        # config.BENCHMARK_LOGGER.info('created cell in: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        

    def _topology(self):
        raise Exception, 'Deprecated'
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'

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
        sim.simdt = 1e-6
        sim.plotdt = 1e-4
        mycell = SupPyrRS(SupPyrRS.prototype, sim.model.path + "/SupPyrRS")
        config.LOGGER.info('Created cell: %s' % (mycell.path))
        vm_table = mycell.soma.insertRecorder('Vm_suppyrrs', 'Vm', sim.data)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=0.4e-9, firstDelay=100e-3, firstWidth=200e-3)
        sim.schedule()
        if mycell.has_cycle():
            config.LOGGER.warning("WARNING!! CYCLE PRESENT IN CICRUIT.")
        t1 = datetime.now()
        sim.run(500e-3)
        t2 = datetime.now()
        delta = t2 - t1
        config.BENCHMARK_LOGGER.info('simulation time: %g' % (delta.seconds + 1e-6 * delta.microseconds))
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            mus_t = linspace(0, sim.simtime*1e3, len(mus_vm))
            try:
                nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_suppyrrs.plot')
                nrn_t = nrn_vm[:, 0]
                nrn_vm = nrn_vm[:, 1]
                # nrn_ca = config.pylab.loadtxt('../nrn/mydata/Ca_suppyrrs.plot')
                # nrn_ca = nrn_ca[:,1]
                config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            except IOError:
                pass
            data = config.pylab.zeros((len(mus_vm), 2))
            data[:, 0] = mus_t[:]
            data[:, 1] = mus_vm[:]
            
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
            # if ca_table:
            #     ca_array = config.pylab.array(ca_table)
            #     config.pylab.plot(nrn_t, -nrn_ca, 'r-', label='nrn (-)ca')
            #     config.pylab.plot(mus_t, -ca_array, 'b-.', label='mus (-)ca')
            #     print config.pylab.amax(ca_table)
            config.pylab.legend()
            config.pylab.title('suppyrRS')
            config.pylab.show()
            data_array = config.pylab.zeros((len(mus_vm), 2))
            data_array[:, 0] = mus_t[:]
            data_array[:, 1] = mus_vm[:]
            config.pylab.savetxt('Vm_suppyrRS.dat', data_array)
        

import unittest
import uuid
class SupPyrRSTestCase(unittest.TestCase):
    def setUp(self):
        self.conductance_densities = defaultdict(list)
        self.conductance_densities['NaF'] = [10.0 * x for x in [400, 187.5, 93.75, 12.5, 12.5, 125, 93.75, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5]]
        self.conductance_densities['NaP'] = [10.0 * x for x in [0, 0.12, 0.06, 0.008, 0.008, 0.08, 0.06, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008]]
        self.conductance_densities['KDR'] = [10.0 * x for x in [400, 125, 93.75, 6.25, 6.25, 125, 93.75, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25]]
        self.conductance_densities['KC'] = [10.0 * x for x in [0, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]]
        self.conductance_densities['KA'] = [10.0 * x for x in [2, 30, 2, 2, 2, 30, 2, 2, 2, 2, 2, 2, 2]]
        self.conductance_densities['KM'] = [10.0 * x for x in [0, 7.5, 7.5,  7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5]]
        self.conductance_densities['K2'] = [10.0 * 0.1] * 13
        self.conductance_densities['KAHP'] = [10.0 * x for x in [0, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]]
        self.conductance_densities['CaL'] = [10.0 * x for x in [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        self.conductance_densities['CaT'] = [10.0 * x for x in [0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]
        self.conductance_densities['AR'] = [10.0 * x for x in [0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]]
        self.sim = Simulation('SupPyrRS')
        path = self.sim.model.path + '/TestSupPyrRS'
        config.LOGGER.debug('Creating cell %s' % path)
        TraubCell.adjust_chanlib(SupPyrRS.chan_params)
        self.cell = SupPyrRS(SupPyrRS.prototype,  "%s/SupPyrRS%d" % (self.sim.model.path, uuid.uuid4().int))
        config.LOGGER.debug('Cell created')
        self.sim.schedule()

    
    def test_compartment_count(self):
        for comp_no in range(SupPyrRS.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))

    def test_initVm(self):
        for comp_no in range(SupPyrRS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -70e-3)

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

    def test_conductances(self):
        for level, comp_nums in self.cell.level.items():
            for comp_num in comp_nums:
                comp = self.cell.comp[comp_num]
                print 'Here'
                for chan_id in moose.context.getWildcardList('%s/#[TYPE=HHChannel]' % (comp.path), True):
                    print chan_id.path()
                    channel = moose.HHChannel(chan_id)
                    channame = channel.name
                    gbar = channel.Gbar / comp.sarea()
                    if level != 0 and comp_num != 1: # compensate for dendritic area doubling for spines
                        gbar /= 2.0 
                    self.assertAlmostEqual(self.conductance_densities[channame][level], gbar)
                
    
# test main --
from simulation import Simulation
from subprocess import call
if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_suppyrRS.hoc'], cwd='../nrn')
    # unittest.main()
    SupPyrRS.test_single_cell()
    

# 
# suppyrrs.py ends here
