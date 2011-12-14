# tcr.py --- 
# 
# Filename: tcr.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Oct 16 10:14:07 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 11:08:03 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 87
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This is a redoing of the Thalamocortical relay cells using prototype file.
# It is a translation of the cell in Traub et al, 2005 model.
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


class TCR(TraubCell):
    chan_params = {
        'ENa': 50e-3,
        'EK': -95e-3,
        'EAR': -35e-3,
        'ECa': 125e-3,
        'EGABA': -81e-3,
        'TauCa': 20e-3,
        'X_AR': 0.25
    }
    ca_dep_chans = ['KAHP_SLOWER', 'KC']
    num_comp = 137
    presyn = 135
    level = TraubCell.readlevels('TCR.levels')
    depth = None
    proto_file = 'TCR.p'
    prototype = TraubCell.read_proto(proto_file, "TCR", level_dict=level, depth_dict=depth, params=chan_params)
    def __init__(self, *args):
        TraubCell.__init__(self, *args)
        moose.CaConc(self.soma.path + '/CaPool').tau = 50e-3
	
    def _topology(self):
        raise Exception, 'Deprecated'
        self.presyn = 135
        self.level[1].add(self.comp[1])
        for ii  in range(2, 120, 13):
            self.level[2].add(self.comp[ii])
        for ii  in range(3, 121, 13):
            self.level[3].add(self.comp[ii])
            self.level[3].add(self.comp[ii+1])
            self.level[3].add(self.comp[ii+2])
        for ii  in range(6, 124, 13):
            for kk in range(0,9):
                self.level[4].add(self.comp[ii+kk])
        for ii  in range(132, 138):
            self.level[0].add(self.comp[ii])
            
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'
        for comp in self.comp[1:]:
            comp.Em = -70e-3
	    comp.initVm = -70e-3

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'
	for comp in self.comp[1:]:
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
#                         if not obj.name in self.chan_list:
#                             obj.Gbar = 0.0
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
	    if ca_pool:
		for channel in ca_chans:
		    channel.connect('IkSrc', ca_pool, 'current')
		    print comp.name, ':', channel.name, 'connected to', ca_pool.name
		for channel in ca_dep_chans:
		    channel.useConcentration = 1
		    ca_pool.connect("concSrc", channel, "concen")
		    print comp.name, ':', ca_pool.name, 'connected to', channel.name


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
        mycell = TCR(TCR.prototype, sim.model.path + "/TCR")
        print 'Created cell:', mycell.path
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_TCR', 'Vm', sim.data)
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
        if config.has_pylab:
            try:
                nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_TCR.plot')
                nrn_t = nrn_vm[:, 0]
                nrn_vm = nrn_vm[:, 1]
                config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            except IOError:
                pass
            mus_vm = config.pylab.array(vm_table) * 1e3
            mus_t = linspace(0, sim.simtime*1e3, len(mus_vm))
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
            config.pylab.legend()
            config.pylab.title('TCR')
            config.pylab.show()
import unittest
import uuid
class TCRTestCase(unittest.TestCase):
    def setUp(self):
        self.conductance_densities = defaultdict(list)
        self.conductance_densities['NaF_TCR'] = [10.0 * x for x in [400, 100, 100, 5, 5]]
        self.conductance_densities['NaPF_TCR'] = [10.0 * x for x in [0.8, 0.2, 0.2, 0.01, 0.01]]
        # The multiplication of g(KDR) by 0.45 for axon, soma and
        # level2 is also implemented in NEURON; mentioned in fortran
        # comment, but not implemented
        self.conductance_densities['KDR'] = [10.0 * x for x in [180, 33.75, 22.5, 0, 0]]
        self.conductance_densities['KC'] = [10.0 * x for x in [0, 12, 12, 20, 20]]
        # multiplication of g(KA) by 0.2 occurs in neuron code, but
        # not mentioned in paper. mentioned in comment pf fortram code
        # also
        self.conductance_densities['KA'] = [10.0 * x for x in [1, 6, 6, 0.2, 0.2]] 
        self.conductance_densities['KM'] = [10.0 * x for x in [0, 0.5, 0.5, 0.5, 0.5]]
        self.conductance_densities['K2'] = [10.0 * x for x in [0.5, 2, 2, 2, 2]]
        self.conductance_densities['KAHP_SLOWER'] = [10.0 * x for x in [0, 0.05, 0.05, 0.05, 0.05]]
        self.conductance_densities['CaL'] = [10.0 * x for x in [0, 0.5, 0.5, 0.25, 0.25]]
        self.conductance_densities['CaT'] = [10.0 * x for x in [0, 0.5, 5, 3, 0.5]]
        self.conductance_densities['AR'] = [10.0 * x for x in [0, 0.25, 0.5, 0.3, 0.3]]
        
        self.sim = Simulation('TCR')
        path = self.sim.model.path + '/TestTCR'
        config.LOGGER.debug('Creating cell %s' % path)
        self.cell = TCR(TCR.prototype,  "%s/TCR%d" % (self.sim.model.path, uuid.uuid4().int))
        print self.cell.comp[1].Em
        config.LOGGER.debug('Cell created')
        self.sim.schedule()

    
    def test_compartment_count(self):
        for comp_no in range(TCR.num_comp):
            path = '%s/comp_%d' % (self.cell.path, comp_no + 1)
            self.assertTrue(config.context.exists(path))

    def test_initVm(self):
        for comp_no in range(TCR.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].initVm, -70e-3)

    def test_Em(self):
        for comp_no in range(TCR.num_comp):
            self.assertAlmostEqual(self.cell.comp[comp_no + 1].Em, -70e-3)

    def test_Ca_connections(self):
        for comp_no in range(TCR.num_comp):
            ca_path = self.cell.comp[comp_no + 1].path + '/CaPool'
            if not config.context.exists(ca_path):
                config.LOGGER.debug('%s : No CaPool' % (self.cell.comp[comp_no + 1].path))
                continue
            caPool = moose.CaConc(ca_path)
            for chan in TCR.ca_dep_chans:
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
        for num in range(TCR.num_comp):
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
                self.assertAlmostEqual(chan.Ek, TCR.chan_params[key])

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
                    # print channel.path
                    self.assertAlmostEqual(self.conductance_densities[channame][level], gbar)
        
        
# test main --
from simulation import Simulation
from subprocess import call
if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_TCR.hoc'], cwd='../nrn')
    TCR.test_single_cell()
    # unittest.main()

# 
# tcr.py ends here
