# supLTS.py --- 
# 
# Filename: supLTS.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Mon Sep 23 00:18:00 2009 (+0530)
# Version: 
# Last-Updated: Mon Dec 12 13:54:03 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 189
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

# Code:
import string
from datetime import datetime
import config
import trbutil
import moose
from cell import *
from capool import CaPool

class SupLTS(TraubCell):
    chan_params = {
        'ENa': 50e-3,
        'EK': -100e-3,
        'ECa': 125e-3,
        'EAR': -40e-3, # dummy
        'EGABA': -75e-3, # Sanchez-Vives et al. 1997 
        'TauCa': 20e-3,
        'X_AR': 0.25
    }
    ca_dep_chans = ['KAHP_SLOWER', 'KC_FAST']
    num_comp = 59
    presyn = 59
    # level maps level number to the set of compartments belonging to it
    level = TraubCell.readlevels("SupLTS.levels")
    # depth stores a map between level number and the depth of the compartments.
    depth = None    
    proto_file = 'SupLTS.p'
    prototype = TraubCell.read_proto(proto_file, "SupLTS", level_dict=level, depth_dict=depth, params=chan_params)
    def __init__(self, *args):
        # start = datetime.now()
        TraubCell.__init__(self, *args)
        caPool = moose.CaConc(self.soma.path + '/CaPool')
        caPool.tau = 50e-3
        # end = datetime.now()
        # delta = end - start
        # config.BENCHMARK_LOGGER.info('created cell in: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        

    def _topology(self):
        raise Exception, 'Deprecated'
	self.presyn = 59

    def _setup_passive(self):
        raise Exception, 'Deprecated'
	for ii in range(SupLTS.num_comp):
            comp = self.comp[1 + ii]
	    comp.initVm = -65e-3

    def _setup_channels(self):
        raise Exception, 'Deprecated'
        for ii in range(SupLTS.num_comp):
            comp = self.comp[1 + ii]
	    ca_pool = None
	    ca_dep_chans = []
	    ca_chans = []
	    for child in comp.children():
		obj = moose.Neutral(child)
		if obj.name == 'CaPool':
		    ca_pool = moose.CaConc(child)
		    ca_pool.tau = 20e-3
		else:
		    obj_class = obj.className
		    if obj_class == 'HHChannel':
			obj = moose.HHChannel(child)
			pyclass = eval(obj.name)
			if issubclass(pyclass, KChannel):
			    obj.Ek = -100e-3
			    if issubclass(pyclass, KCaChannel):
				ca_dep_chans.append(obj)
			elif issubclass(pyclass, NaChannel):
			    obj.Ek = 50e-3
			elif issubclass(pyclass, CaChannel):
			    obj.Ek = 125e-3
			    if issubclass(pyclass, CaL):
				ca_chans.append(obj)
			elif issubclass(pyclass, AR):
			    obj.Ek = -40e-3
	    if ca_pool:
		for channel in ca_chans:
		    channel.connect('IkSrc', ca_pool, 'current')
		    print comp.name, ':', channel.name, 'connected to', ca_pool.name
		for channel in ca_dep_chans:
		    channel.useConcentration = 1
		    ca_pool.connect("concSrc", channel, "concen")
		    print comp.name, ':', ca_pool.name, 'connected to', channel.name

	obj = moose.CaConc(self.soma.path + '/CaPool')
        obj.tau = 50e-3
        
    def check_model(self):
        self.check_compartment_count()
        self.check_initVm()
        self.check_Em()
        self.check_Ca_connection()
        self.check_reversal_potentials()

    def check_compartment_count(self):
        for comp_no in range(SupLTS.num_comp):
            path = '%s/comp_%d' % (self.path, comp_no + 1)
            assert (config.context.exists(path))

    def check_initVm(self):
        for comp_no in range(SupLTS.num_comp):
            assert trbutil.almost_equal(self.cell.comp[comp_no + 1].initVm, -65e-3)
    
    def check_Em(self):
        for comp_no in range(SupLTS.num_comp):
            assert trbutil.almost_equal(self.cell.comp[comp_no + 1].Em, -65e-3)
        
    def check_Ca_connection(self):
        for comp_no in range(SupLTS.num_comp):
            ca_path = self.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                continue
            caPool = moose.CaConc(ca_path)
            for chan in SupLTS.ca_dep_chans:
                chan_path = self.comp[comp_no + 1].path + '/' + chan
                if not config.context.exists(chan_path):
                    continue
                chan_obj = moose.HHChannel(chan_path)
                assert (len(chan_obj.neighbours('concen')) > 0)
            sources = caPool.neighbours('current')
            assert (len(sources) > 0)
            for chan in sources:
                assert (chan.path().endswith('CaL'))
        
    def check_reversal_potentials(self):
        for num in range(SupLTS.num_comp):
            comp = self.comp[num + 1]
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
                assert trbutil.almost_equal(chan.Ek, SupLTS.chan_params[key])
        
        

    @classmethod
    def test_single_cell(cls):
        """Simulates a single superficial LTS cell and plots the Vm and [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = SupLTS(SupLTS.prototype, sim.model.path + "/SupLTS")
        print 'Created cell:', mycell.path
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_supLTS', 'Vm', sim.data)
        ca_conc_path = mycell.soma.path + '/CaPool'
        ca_table = None
        if config.context.exists(ca_conc_path):
            ca_conc = moose.CaConc(ca_conc_path)
            ca_table = moose.Table('Ca_supLTS', sim.data)
            ca_table.stepMode = 3
            ca_conc.connect('Ca', ca_table, 'inputRequest')
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=50e-3, firstWidth=50e-3)
        sim.schedule()
        if mycell.has_cycle():
            print "WARNING!! CYCLE PRESENT IN CICRUIT."
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        print 'simulation time: ', delta.seconds + 1e-6 * delta.microseconds
        sim.dump_data('data')
        mycell.dump_cell('supLTS.txt')
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_supLTS.plot')
            nrn_t = nrn_vm[:, 0]
            mus_t = linspace(0, sim.simtime*1e3, len(mus_vm))
            nrn_vm = nrn_vm[:, 1]
            nrn_ca = config.pylab.loadtxt('../nrn/mydata/Ca_supLTS.plot')
            nrn_ca = nrn_ca[:,1]
            title = 'SupLTS:'
            config.pylab.title(title)
            config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
            config.pylab.legend()
            config.pylab.show()


import unittest
import uuid

class SupLTSTestCase(unittest.TestCase):
    def setUp(self):
        self.conductance_densities = defaultdict(list)
        self.conductance_densities['NaF2'] = [10.0 * x for x in [400, 60, 60, 60, 10, 10, 10, 10, 10, 10]]
        
        self.conductance_densities['KDR_FS'] = [10.0 * x for x in [400, 100, 100, 100, 10, 10, 10, 10, 10, 10]]
        self.conductance_densities['KC_FAST'] = [10.0 * x for x in [0, 25, 25, 25, 25, 25, 25, 25, 25, 25]]
        self.conductance_densities['KA'] = [10.0 * 1] * 10
        self.conductance_densities['KM'] = [10 * x for x in [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]
        self.conductance_densities['K2'] = [10.0 * 0.5] * 10
        self.conductance_densities['KAHP_SLOWER'] = [10.0 * x for x in [0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]
        self.conductance_densities['CaT_A'] = [10.0 * x for x in [0.0, 0.05, 0.05, 0.05, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]]
        self.conductance_densities['CaL'] = [10.0 * x for x in [0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]]

        self.sim = Simulation('SupLTS')
        path = self.sim.model.path + '/' + 'TestSupLTS'
        config.LOGGER.debug('Creating cell %s' % (path))
        TraubCell.adjust_chanlib(SupLTS.chan_params)
        self.cell = SupLTS(SupLTS.prototype, '%s/SupLTS%d' % (self.sim.model.path, uuid.uuid4().int))
        config.LOGGER.debug('Cell created')
        for handler in config.LOGGER.handlers:
            handler.flush()
        self.sim.schedule()

    def test_compartment_count(self):
        for comp_no in range(SupLTS.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))
    
    def test_initVm(self):
        for comp_no in range(SupLTS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -65e-3)

    def test_Em(self):
        for comp_no in range(SupLTS.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].Em, -65e-3)

    def test_Ca_connections(self):
        for comp_no in range(SupLTS.num_comp):
            ca_path = self.cell.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                continue
            caPool = moose.CaConc(ca_path)
            for chan in SupLTS.ca_dep_chans:
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
        for num in range(SupLTS.num_comp):
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
                self.assertAlmostEqual(chan.Ek, SupLTS.chan_params[key])
                    
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
                self.assertAlmostEqual(self.conductance_densities[channame][level], gbar)
                        
        
# test main --
from simulation import Simulation
from subprocess import call

if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_supLTS.hoc'], cwd='../nrn')
    SupLTS.test_single_cell()
    # unittest.main()



# 
# supLTS.py ends here
