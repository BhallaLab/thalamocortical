# spinystellate.py --- 
# 
# Filename: spinystellate.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Tue Sep 29 11:43:22 2009 (+0530)
# Version: 
# Last-Updated: Thu Oct 27 14:15:45 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 530
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
        'EAR': -40e-3,
        'ECa': 100e-3,
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
        # for arg in args: print arg
        TraubCell.__init__(self, *args)
        # print 'TraubCell.__init__ passed'
        # print 'Soma:', self.soma.path
        soma_ca_pool = moose.CaConc(self.soma.path + '/CaPool')
        # print 'CaPool generated'
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


        mycell.soma.x0 = 0.0
        mycell.soma.y0 = 0.0
        mycell.soma.z0 = 0.0
        mycell.soma.x = 0.0
        mycell.soma.y = 0.0
        mycell.soma.z = mycell.soma.length
        # mycellview = MyCellView(mycell)
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
            mus_t = linspace(0, sim.simtime, len(mus_vm))
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
            config.pylab.show()

        

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
