# nRT.py --- 
# 
# Filename: nRT.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Oct 16 15:18:24 2009 (+0530)
# Version: 
# Last-Updated: Wed Feb 17 17:22:09 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 53
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
import config
import trbutil
import moose
from cell import *
from capool import CaPool


class nRT(TraubCell):
    prototype = TraubCell.read_proto("nRT.p", "nRT")
    ca_dep_chans = ['KAHP','KAHP_SLOWER', 'KAHP_DP', 'KC', 'KC_FAST']
    unblocked_chans = [ 'CaL', 'KAHP_SLOWER']#['NaPF', 'KDR_FS', 'NaF2_nRT', 'K2', 'KM']
    def __init__(self, *args):
	TraubCell.__init__(self, *args)
	
    def _topology(self):
        self.presyn = 59
        self.level[1].add(self.comp[1])
        for ii in range(2, 42, 13):
            self.level[2].add(self.comp[ii])
        for ii in range(3, 43, 13):
            self.level[3].add(self.comp[ii])
            self.level[3].add(self.comp[ii + 1])
        for ii in range(5, 45, 13):
            self.level[4].add(self.comp[ii])
            self.level[4].add(self.comp[ii + 1])
            self.level[4].add(self.comp[ii + 2])
        for ii in range(8, 47, 13):
            self.level[5].add(self.comp[ii])
            self.level[5].add(self.comp[ii + 1])
            self.level[5].add(self.comp[ii + 2])
        for ii in range(11, 51, 13):
            self.level[6].add(self.comp[ii])
            self.level[7].add(self.comp[ii + 1])
            self.level[8].add(self.comp[ii + 2])
            self.level[9].add(self.comp[ii + 3])
        for ii in range(54, 60):
            self.level[0].add(self.comp[ii])
    
    def _setup_passive(self):
        for comp in self.comp[1:]:
	    comp.initVm = -75e-3

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
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
                        gbar = obj.Gbar
#                         if not obj.name in nRT.unblocked_chans:
#                             obj.Gbar = 0.0
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
                            obj.X = 0.0
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


    @classmethod
    def test_single_cell(cls):
        """Simulates a single nRT cell and plots the Vm and [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = nRT(nRT.prototype, sim.model.path + "/nRT")
        config.LOGGER.info('Created cell: %s' % (mycell.path))
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_nRT', 'Vm', sim.data)
        ca_conc_path = mycell.soma.path + '/CaPool'
        ca_table = None
        if config.context.exists(ca_conc_path):
            ca_conc = moose.CaConc(ca_conc_path)
            ca_table = moose.Table('Ca_nRT', sim.data)
            ca_table.stepMode = 3
            ca_conc.connect('Ca', ca_table, 'inputRequest')
        kc_path = mycell.soma.path + '/KC'
        gk_table = None
        if config.context.exists(kc_path):
            gk_table = moose.Table('gkc', sim.data)
            gk_table.stepMode = 3
            kc = moose.HHChannel(kc_path)
            kc.connect('Gk', gk_table, 'inputRequest')
            pymoose.showmsg(ca_conc)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=50.0e-3, firstWidth=50e-3)
#         pulsegen1 = mycell.soma.insertPulseGen('pulsegen1', sim.model, firstLevel=3e-7, firstDelay=150e-3, firstWidth=10e-3)

        sim.schedule()
        if mycell.has_cycle():
            config.LOGGING.warning("WARNING!! CYCLE PRESENT IN CICRUIT.")
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        config.LOGGER.info('simulation time: %g'  % (delta.seconds + 1e-6 * delta.microseconds))
        sim.dump_data('data')
        mycell.dump_cell('nRT.txt')
        
        mus_vm = pylab.array(vm_table) * 1e3
        nrn_vm = pylab.loadtxt('../nrn/mydata/Vm_nRT.plot')
        nrn_t = nrn_vm[:, 0]
        mus_t = linspace(0, nrn_t[-1], len(mus_vm))
        nrn_vm = nrn_vm[:, 1]
        nrn_ca = pylab.loadtxt('../nrn/mydata/Ca_nRT.plot')
        nrn_ca = nrn_ca[:,1]
        pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
        pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
#         if ca_table:
#             ca_array = pylab.array(ca_table)
#             pylab.plot(nrn_t, -nrn_ca, 'r-', label='nrn (-)ca')
#             pylab.plot(mus_t, -ca_array, 'b-.', label='mus (-)ca')
#             print pylab.amax(ca_table)
        pylab.legend()
        pylab.show()
        
        
# test main --
from simulation import Simulation
import pylab
from subprocess import call
if __name__ == "__main__":
    call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_nRT.hoc'], cwd='../nrn')
    nRT.test_single_cell()




# 
# nRT.py ends here
