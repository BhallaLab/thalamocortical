# deepLTS.py --- 
# 
# Filename: deepLTS.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Oct 16 19:32:34 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 11:04:45 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 52
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


class DeepLTS(TraubCell):
    chan_params = {
        'ENa': 50e-3,
        'EK': -100e-3,
        'EAR': -40e-3,
        'ECa': 125e-3,
        'EGABA': -75e-3, # Sanchez-Vives et al. 1997 
        'TauCa': 20e-3,
        'X_AR': 0.25
    }
    ca_dep_chans = ['KAHP_SLOWER', 'KC_FAST']
    num_comp = 59
    presyn = 59
    # level maps level number to the set of compartments belonging to it
    level = TraubCell.readlevels("DeepLTS.levels")
    # depth stores a map between level number and the depth of the compartments.
    depth = None    

    proto_file = 'DeepLTS.p'
    prototype = TraubCell.read_proto(proto_file, "DeepLTS", level_dict=level, depth_dict=depth, params=chan_params)
    def __init__(self, *args):
        TraubCell.__init__(self, *args)
        moose.CaConc(self.soma.path + '/CaPool').tau = 50e-3
	
    def _topology(self):
        raise Exception, 'Deprecated'
        self.presyn = 59
        self.level[1].add(self.comp[1])
        for ii in range(2, 42, 13):
            self.level[2].add(self.comp[ii])
        for ii in range(3, 43, 13):
            self.level[3].add(self.comp[ii])
            self.level[3].add(self.comp[ii+1])
        for ii in range(5, 45, 13):
            self.level[4].add(self.comp[ii])
            self.level[4].add(self.comp[ii+1])
            self.level[4].add(self.comp[ii+2])
        for ii in range(8, 48, 13):
            self.level[5].add(self.comp[ii])
            self.level[5].add(self.comp[ii+1])
            self.level[5].add(self.comp[ii+2])
        for ii in range(11, 51, 13):
            self.level[6].add(self.comp[ii])
            self.level[7].add(self.comp[ii+1])
            self.level[8].add(self.comp[ii+2])
            self.level[9].add(self.comp[ii+3])
        for ii in range(54, 60):
            self.level[0].add(self.comp[ii])
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'


    @classmethod
    def test_single_cell(cls):
        """Simulates a single deep LTS cell and plots the Vm and [Ca2+]"""
        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = DeepLTS(DeepLTS.prototype, sim.model.path + "/DeepLTS")
        
        config.LOGGER.debug(('Created cell: %s' % mycell.path))
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_deepLTS', 'Vm', sim.data)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=50e-3, firstWidth=50e-3)

        sim.schedule()
        if mycell.has_cycle():
            config.LOGGER.warning("WARNING!! CYCLE PRESENT IN CICRUIT.")

        sim.run(200e-3)
        sim.dump_data('data')
        mycell.dump_cell('deepLTS.txt')
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            mus_t = linspace(0, sim.simtime * 1e3, len(mus_vm))
            try:
                nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_supLTS.plot')
                nrn_t = nrn_vm[:, 0]
                nrn_vm = nrn_vm[:, 1]
                config.pylab.plot(nrn_t, nrn_vm, 'y-', label='nrn vm')
            except IOError:
                print 'NEURON Data not available.'
            config.pylab.plot(mus_t, mus_vm, 'g-.', label='mus vm')
            config.pylab.legend()
            config.pylab.title('DeepLTS')
            config.pylab.show()

        
# test main --
from simulation import Simulation
from subprocess import call
if __name__ == "__main__":
    DeepLTS.test_single_cell()





# 
# deepLTS.py ends here
