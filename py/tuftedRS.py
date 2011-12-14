# tuftedRS.py --- 
# 
# Filename: tuftedRS.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Oct 16 13:42:14 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 11:08:54 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 68
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

from datetime import datetime
import moose
import config
import trbutil
from cell import *
from capool import CaPool


class TuftedRS(TraubCell):
    chan_params = {
        'ENa': 50e-3,
        'EK': -95e-3,
        'EAR': -35e-3,
        'ECa': 125e-3,
        'EGABA': -75e-3, # Sanchez-Vives et al. 1997 
        'TauCa': 1e-3/0.075,
        'X_AR': 0.25
        }
    ca_dep_chans = ['KAHP_DP', 'KC']
    num_comp = 61
    presyn = 60
    level = TraubCell.readlevels('TuftedRS.levels')
    depth = {
        1: 1800 * 1e-6,
        2: 1845 * 1e-6,
        3: 1890 * 1e-6,
        4: 1935 * 1e-6,
        5: 1760 * 1e-6,
        6: 1685 * 1e-6,
        7: 1610 * 1e-6,
        8: 1535 * 1e-6,
        9: 1460 * 1e-6,
        10: 1385 * 1e-6,
        11: 1310 * 1e-6,
        12: 1235 * 1e-6,
        13: 1160 * 1e-6,
        14: 1085 * 1e-6,
        15: 1010 * 1e-6,
        16: 935 * 1e-6,
        17: 860 * 1e-6,
        18: 790 * 1e-6,
        }
    proto_file = 'TuftedRS.p'
    prototype = TraubCell.read_proto(proto_file, "TuftedRS", level_dict=level, depth_dict=depth, params=chan_params)
    ca_dep_chans = ['KAHP_DP', 'KC']

    def __init__(self, *args):
        TraubCell.__init__(self, *args)
        moose.CaConc(self.soma.path + '/CaPool').tau = 100e-3
        # Special case: individually specified beta_cad's in level  2
        moose.CaConc(self.comp[2].path + '/CaPool').tau  =   1e-3/0.02
        moose.CaConc(self.comp[5].path + '/CaPool' ).tau = 1e-3 /  0.02
        moose.CaConc(self.comp[6].path + '/CaPool' ).tau = 1e-3 /  0.02
	
    def _topology(self):
        raise Exception, 'Deprecated'
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'


    @classmethod
    def test_single_cell(cls):
        """Simulates a single thalamocortical relay cell
        and plots the Vm and [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = TuftedRS(TuftedRS.prototype, sim.model.path + "/TuftedRS")
        print 'MOOSE: Created cell:', mycell.path
        vm_table = mycell.comp[cls.presyn].insertRecorder('Vm_tuftRS', 'Vm', sim.data)
        # ca_conc_path = mycell.soma.path + '/CaPool'
        # ca_table = None
        # if config.context.exists(ca_conc_path):
        #     ca_conc = moose.CaConc(ca_conc_path)
        #     ca_table = moose.Table('Ca_tuftRS', sim.data)
        #     ca_table.stepMode = 3
        #     ca_conc.connect('Ca', ca_table, 'inputRequest')
        # kc_path = mycell.soma.path + '/KC'
        # gk_table = None
        # if config.context.exists(kc_path):
        #     gk_table = moose.Table('gkc', sim.data)
        #     gk_table.stepMode = 3
        #     kc = moose.HHChannel(kc_path)
        #     kc.connect('Gk', gk_table, 'inputRequest')
        #     pymoose.showmsg(ca_conc)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=10e-10, firstDelay=0.0, firstWidth=50e-3)
#         pulsegen1 = mycell.soma.insertPulseGen('pulsegen1', sim.model, firstLevel=3e-7, firstDelay=150e-3, firstWidth=10e-3)

        sim.schedule()
        if mycell.has_cycle():
            print "WARNING!! CYCLE PRESENT IN CICRUIT."
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        print 'MOOSE: simulation time: ', delta.seconds + 1e-6 * delta.microseconds
        sim.dump_data('data')
        # mycell.dump_cell('tuftRS.txt')
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_tuftRS.plot')
            nrn_t = nrn_vm[:, 0]
            mus_t = linspace(0, nrn_t[-1], len(mus_vm))
            nrn_vm = nrn_vm[:, 1]
            # nrn_ca = config.pylab.loadtxt('../nrn/mydata/Ca_tuftRS.plot')
            # nrn_ca = nrn_ca[:,1]
            config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
    #         if ca_table:
    #             ca_array = config.pylab.array(ca_table)
    #             config.pylab.plot(nrn_t, -nrn_ca, 'r-', label='nrn (-)ca')
    #             config.pylab.plot(mus_t, -ca_array, 'b-.', label='mus (-)ca')
    #             print config.pylab.amax(ca_table)
            config.pylab.legend()
            config.pylab.title('tuftedRS')
            config.pylab.show()
        
import unittest
import uuid
class TuftedRSTestCase(unittest.TestCase):
    def setUp(self):
        self.conductance_densities = defaultdict(list)
        self.conductance_densities['NaF'] = [10.0 * x for x in [450, 200, 75, 15, 15, 150, 75, 15, 15, 15, 15, 15, 15, 15, 15, 3, 3, 3, 3, 3]]
        self.conductance_densities['NaP'] = [10.0 * x for x in [0, 0.16, 0.06, 0.012, 0.012, 0.12, 0.06, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.0024, 0.0024, 0.0024, 0.0024, 0.0024]]
        self.conductance_densities['KDR'] = [10.0 * x for x in [450, 170, 75, 0, 0, 120, 75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        self.conductance_densities['KC'] = [10.0 * x for x in [0, 28.8, 28.8, 0.9, 0.9, 28.8, 28.8, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 2.16, 2.16, 2.16, 2.16, 2.16]]
        self.conductance_densities['KA'] = [10.0 * x for x in [0.6, 20, 8, 0.6, 0.6, 8, 8, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]]
        self.conductance_densities['KM'] = [10.0 * x for x in [30, 8.5, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 4, 4, 4, 4, 4]] 
        self.conductance_densities['K2'] = [10.0 * 0.5] * 19
        self.conductance_densities['KAHP_DP'] = [10.0 * x for x in [0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]]
        self.conductance_densities['CaL'] = [10.0 * x for x in [0, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.08, 0.24]]
        self.conductance_densities['CaT'] = [10.0 * x for x in [0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]
        self.conductance_densities['AR'] = [10.0 * x for x in [0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2]]
        self.sim = Simulation('TuftedRS')
        path = self.sim.model.path + '/TestTuftedRS'
        config.LOGGER.debug('Creating cell %s' % path)
        TraubCell.adjust_chanlib(TuftedRS.chan_params)
        self.cell = TuftedRS(TuftedRS.prototype,  "%s/TuftedRS%d" % (self.sim.model.path, uuid.uuid4().int))
        config.LOGGER.debug('Cell created')
        self.sim.schedule()

    
    def test_compartment_count(self):
        for comp_no in range(TuftedRS.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))

    def test_initVm(self):
        for comp_no in range(TuftedRS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -70e-3)

    def test_Em(self):
        for comp_no in range(TuftedRS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].Em, -70e-3)

    def test_Ca_connections(self):
        for comp_no in range(TuftedRS.num_comp):
            ca_path = self.cell.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                config.LOGGER.debug('%s : No CaPool' % (self.cell.comp[comp_no + 1].path))
                continue
            caPool = moose.CaConc(ca_path)
            for chan in TuftedRS.ca_dep_chans:
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
        for num in range(TuftedRS.num_comp):
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
                self.assertAlmostEqual(chan.Ek, TuftedRS.chan_params[key])

    def test_conductances(self):
        for level, comp_nums in self.cell.level.items():
            for comp_num in comp_nums:
                comp = self.cell.comp[comp_num]
                print 'Here'
                for chan_id in moose.context.getWildcardList('%s/#[TYPE=HHChannel]' % (comp.path), True):
                    channel = moose.HHChannel(chan_id)
                    channame = channel.name
                    print channel.path, level
                    gbar = channel.Gbar / comp.sarea()
                    if level != 0 and comp_num != 1: # compensate for dendritic area doubling for spines
                        gbar /= 2.0
                    if channame == 'CaL' and level == 18 and comp_num > 49:
                        self.assertAlmostEqual(self.conductance_densities[channame][19], gbar)
                    else:
                        self.assertAlmostEqual(self.conductance_densities[channame][level], gbar)
        
# test main --
from simulation import Simulation
from subprocess import call
if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_tuftRS.hoc'], cwd='../nrn')
    TuftedRS.test_single_cell()
    # unittest.main()



# 
# tuftedRS.py ends here
