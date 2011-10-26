# tuftedIB.py --- 
# 
# Filename: tuftedIB.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Oct 16 11:44:48 2009 (+0530)
# Version: 
# Last-Updated: Wed Oct 26 11:54:16 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 125
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


class TuftedIB(TraubCell):
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
    level = TraubCell.readlevels('TuftedIB.levels')
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
        18: 790 * 1e-6}
    
    proto_file = 'TuftedIB.p'
    prototype = TraubCell.read_proto(proto_file, "TuftedIB", level_dict=level, depth_dict=depth, params=chan_params)
    def __init__(self, *args):
        TraubCell.__init__(self, *args)
        soma_ca_pool = moose.CaConc(self.soma.path + '/CaPool')
        soma_ca_pool.tau = 100e-3
	# Special case: individually specified beta_cad's in level  2
	moose.CaConc(self.comp[2].path + '/CaPool').tau  =   1e-3/0.02
        moose.CaConc(self.comp[ 5].path + '/CaPool' ).tau = 1e-3 /  0.02
        moose.CaConc(self.comp[ 6].path + '/CaPool' ).tau = 1e-3 /  0.02
	
    def _topology(self):
        raise Exception, 'Deprecated'
    
    def _setup_passive(self):
        raise Exception, 'Deprecated'

    def _setup_channels(self):
        """Set up connections between compartment and channels, and Ca pool"""
        raise Exception, 'Deprecated'


    @classmethod
    def test_single_cell(cls):
        """Simulates a single tufted intrinsically bursting cell and
        plots the Vm and [Ca2+]"""

        config.LOGGER.info("/**************************************************************************")
        config.LOGGER.info(" *")
        config.LOGGER.info(" * Simulating a single cell: %s" % (cls.__name__))
        config.LOGGER.info(" *")
        config.LOGGER.info(" **************************************************************************/")
        sim = Simulation(cls.__name__)
        mycell = TuftedIB(TuftedIB.prototype, sim.model.path + "/TuftedIB")
        print 'Created cell:', mycell.path
        mycell.scale_conductance('CaL', 0.0)
        mycell.scale_conductance('CaT', 0.0)
        mycell.scale_conductance('K2', 0.0)
        # mycell.scale_conductance('KA_IB', 0.0)
        mycell.scale_conductance('KAHP_DP', 0.0)
        mycell.scale_conductance('KC', 0.0)
        mycell.scale_conductance('KDR', 0.0)
        mycell.scale_conductance('KM', 0.0)
        mycell.scale_conductance('NaF', 0.0)
        mycell.scale_conductance('NaP', 0.0)
        mycell.scale_conductance('AR', 0.0)
        for chanid in moose.context.getWildcardList(mycell.path+'/##[TYPE=HHChannel]', True):
            chan = moose.HHChannel(chanid)
            if chan.Gbar != 0.0:
                print chan.path, 'Nonzero Gbar'
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_tuftIB', 'Vm', sim.data)
        ca_conc_path = mycell.soma.path + '/CaPool'
        ca_table = None
        if config.context.exists(ca_conc_path):
            ca_conc = moose.CaConc(ca_conc_path)
            ca_table = moose.Table('Ca_tuftIB', sim.data)
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
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=-0.3e-9, firstDelay=50e-3, firstWidth=20e-3, secondLevel=0.9e-9, secondDelay=100e-3, secondWidth=20e-3)
#         pulsegen1 = mycell.soma.insertPulseGen('pulsegen1', sim.model, firstLevel=3e-7, firstDelay=150e-3, firstWidth=10e-3)
        pulse_table = moose.Table('/data/injection')
        pulse_table.stepMode = 3
        pulse_table.connect('inputRequest', pulsegen, 'output')

        sim.schedule()
        if mycell.has_cycle():
            print "WARNING!! CYCLE PRESENT IN CICRUIT."
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        print 'simulation time: ', delta.seconds + 1e-6 * delta.microseconds
        sim.dump_data('data')
        mycell.dump_cell('tuftIB.txt')
        if config.has_pylab:
            mus_vm = config.pylab.array(vm_table) * 1e3
            nrn_vm = config.pylab.loadtxt('../nrn/mydata/Vm_tuftIB.plot')
            nrn_t = nrn_vm[:, 0]
            mus_t = linspace(0, nrn_t[-1], len(mus_vm))
            nrn_vm = nrn_vm[:, 1]
            nrn_ca = config.pylab.loadtxt('../nrn/mydata/Ca_tuftIB.plot')
            nrn_ca = nrn_ca[:,1]
            config.pylab.plot(nrn_t, nrn_vm, 'r-', label='nrn vm')
            config.pylab.plot(mus_t, mus_vm, 'k--', label='mus vm')
            config.pylab.plot(mus_t, config.pylab.array(pulse_table)*1e10-100, 'b-', label='I_inject(x10 nA)')
    #         if ca_table:
    #             ca_array = config.pylab.array(ca_table)
    #             config.pylab.plot(nrn_t, -nrn_ca, 'r-', label='nrn (-)ca')
    #             config.pylab.plot(mus_t, -ca_array, 'b-.', label='mus (-)ca')
    #             print config.pylab.amax(ca_table)
            config.pylab.legend(loc=0)
            config.pylab.show()
        
        
# test main --
from simulation import Simulation
from subprocess import call
if __name__ == "__main__":
    # call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_tuftIB.hoc'], cwd='../nrn')
    TuftedIB.test_single_cell()



# 
# tuftedIB.py ends here
