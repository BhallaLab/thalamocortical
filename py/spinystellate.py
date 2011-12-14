# spinystellate.py --- 
# 
# Filename: spinystellate.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Tue Sep 29 11:43:22 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 11:05:53 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 580
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
# Code:

from datetime import datetime
import config
import trbutil
import moose
from cell import *
from capool import CaPool
# from cellview import MyCellView

class SpinyStellate(TraubCell):
    chan_params = {
        'ENa': 50e-3,
        'EK': -100e-3,
        'EAR': -40e-3, # increasing EAR brings the spikes earlier, -37.5 gives an exact match with NEURON model.
        'ECa': 125e-3,
        'EGABA': -75e-3, # Sanchez-Vives et al. 1997 
        'TauCa': 20e-3,
        'X_AR': 0.0
        }
    ca_dep_chans = ['KAHP_SLOWER', 'KC_FAST']
    num_comp = 59
    presyn = 57
    level = TraubCell.readlevels('SpinyStellate.levels')
    depth = None
    proto_file = 'SpinyStellate.p'
    prototype = TraubCell.read_proto("SpinyStellate.p", "SpinyStellate", level_dict=level, depth_dict=depth, params=chan_params)
    def __init__(self, *args):
        # start = datetime.now()
        TraubCell.__init__(self, *args)
        soma_ca_pool = moose.CaConc(self.soma.path + '/CaPool')
        soma_ca_pool.tau = 50e-3
        # end = datetime.now()
        # delta = end - start
        # config.BENCHMARK_LOGGER.info('created cell in: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))

    def _topology(self):
        raise Exception, 'Deprecated method.'

    def _setup_passive(self):
        raise Exception, 'Deprecated. All passive properties including initVm and Em are set in .p file.'
    

    def _setup_channels(self):
        """Set up connection between CaPool, Ca channels, Ca dependnet channels."""
        raise Exception, 'Deprecated.'

    @classmethod
    def test_single_cell(cls):
        """Simulates a single spiny stellate cell and plots the Vm and
        [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        proto = TraubCell.read_proto("SpinyStellate.p", "SpinyStellate_1", SpinyStellate.chan_params)
        print 'Model path', sim.model.path , 'proto', proto

        mycell = SpinyStellate(proto, sim.model.path + "/SpinyStellate")
        print 'Cell created'
        for handler in config.LOGGER.handlers:
            handler.flush()


        config.LOGGER.debug('Created cell: %s' % (mycell.path))
        # for neighbour in mycell.soma.neighbours('raxial'):
        #     print 'RAXIAL', neighbours.path()
        # for neighbour in mycell.soma.neighbours('axial'):
        #     print 'AXIAL', neighbour.path()
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_spinstell', 'Vm', sim.data)
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
        for msg in moose.Neutral('/model/SpinyStellate/solve').inMessages():
            print msg
        for msg in moose.Neutral('/model/SpinyStellate/solve').outMessages():
            print msg
        # sim.dump_data('data')
        # mycell.dump_cell('spinstell.txt')
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            mus_t = linspace(0, sim.simtime*1e3, len(mus_vm))
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
            try:
                nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_spinstell.plot')
                nrn_t = nrn_vm[:, 0]
                nrn_vm = nrn_vm[:, 1]
                nrn_ca = config.pylab.loadtxt('../nrn/mydata/Ca_spinstell.plot')
                config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            except IOError:
                pass
            config.pylab.legend(loc=0)
            config.pylab.title('spinystellate')
            config.pylab.show()

        
import unittest
import uuid
class SpinyStellateTestCase(unittest.TestCase):
    def setUp(self):
        self.conductance_densities = defaultdict(list)
        self.conductance_densities['NaF2'] = [10.0 * x for x in [400, 150, 75, 75, 5, 5, 5, 5, 5, 5]]
        # Strange: the paper says 400 mS/cm^2 density of NaP in level
        # 0 in table A3, but in NEURON code it does not exist for
        # level 0.
        self.conductance_densities['NaPF_SS'] = [10.0 * x for x in [0.4, 0.15, 0.075, 0.075, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]]
        self.conductance_densities['KDR_FS'] = [10.0 * x for x in [400, 100, 75, 75, 0, 0, 0, 0, 0, 0]]
        self.conductance_densities['KC_FAST'] = [10.0 * x for x in [0, 10, 10, 10, 10, 0, 0, 0, 0, 0]]
        self.conductance_densities['KA'] = [10.0 * x for x in [2, 30, 30, 2, 2, 2, 2, 2, 2, 2]]
        self.conductance_densities['KM'] = [10.0 * x for x in [0, 3.75, 3.75,  3.75, 3.75, 3.75, 3.75, 3.75, 3.75, 3.75]]
        self.conductance_densities['K2'] = [10.0 * 0.1] * 10
        self.conductance_densities['KAHP_SLOWER'] = [10.0 * x for x in [0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]
        self.conductance_densities['CaL'] = [10.0 * x for x in [0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 3, 3, 3]]
        self.conductance_densities['CaT_A'] = [10.0 * x for x in [0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]
        self.conductance_densities['AR'] = [10.0 * x for x in [0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]]
        self.sim = Simulation('SpinyStellate')
        path = self.sim.model.path + '/TestSpinyStellate'
        config.LOGGER.debug('Creating cell %s' % path)
        TraubCell.adjust_chanlib(SpinyStellate.chan_params)
        self.cell = SpinyStellate(SpinyStellate.prototype,  "%s/SpinyStellate%d" % (self.sim.model.path, uuid.uuid4().int))
        config.LOGGER.debug('Cell created')
        self.sim.schedule()

    
    def test_compartment_count(self):
        for comp_no in range(SpinyStellate.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))

    def test_initVm(self):
        for comp_no in range(SpinyStellate.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -65e-3)

    def test_Em(self):
        for comp_no in range(SpinyStellate.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].Em, -65e-3)

    def test_Ca_connections(self):
        for comp_no in range(SpinyStellate.num_comp):
            ca_path = self.cell.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                config.LOGGER.debug('%s : No CaPool' % (self.cell.comp[comp_no + 1].path))
                continue
            caPool = moose.CaConc(ca_path)
            for chan in SpinyStellate.ca_dep_chans:
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
        for num in range(SpinyStellate.num_comp):
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
                self.assertAlmostEqual(chan.Ek, SpinyStellate.chan_params[key])

    def test_conductances(self):
        for level, comp_nums in self.cell.level.items():
            for comp_num in comp_nums:
                comp = self.cell.comp[comp_num]
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
#     call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_spinstell.hoc'], cwd='../nrn')
    SpinyStellate.test_single_cell()
    # unittest.main()
    # test_creation_time(1000)

# 
# spinystellate.py ends here
