# suppyrrs.py --- 
# 
# Filename: suppyrrs.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Aug  7 13:59:30 2009 (+0530)
# Version: 
# Last-Updated: Mon Aug 31 21:43:28 2009 (+0530)
#           By: subhasis ray
#     Update #: 201
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
import moose
import config
from cell import *
from capool import CaPool

class SupPyrRS(TraubCell):
    prototype = TraubCell.read_proto("SupPyrRS.p", "SupPyrRS")
    def __init__(self, *args):
	TraubCell.__init__(self, *args)
	
    def _topology(self):
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
        pass

    def _setup_channels(self):
        for i in range(len(self.level)):
            for comp in self.level[i]:
                for child in comp.children():
                    obj = moose.Neutral(child)
                    if obj.name == 'CaPool':
                        obj = moose.CaConc(child)
                        obj.tau = 1e-3/0.05
                        print obj.path, 'set tau to', obj.tau
                    else:
                        obj_class = obj.className
                        if obj_class == "HHChannel":
                            obj = moose.HHChannel(child)
                            pyclass = eval(obj.name)
                            if issubclass(pyclass, KChannel):
                                obj.Ek = -95e-3
                            elif issubclass(pyclass, NaChannel):
                                obj.Ek = 50e-3
                            elif issubclass(pyclass, CaChannel):
                                obj.Ek = 125e-3
                            elif issubclass(pyclass, AR):
                                obj.Ek = -35e-3
        obj = moose.CaConc(self.soma.path + '/CaPool')
        obj.tau = 1e-3 / 0.01
        print obj.path, 'set tau to', obj.tau

    @classmethod
    def test_single_cell(cls):
        sim = Simulation()
        mycell = SupPyrRS(SupPyrRS.prototype, sim.model.path + "/SupPyrRS")
        print 'Created cell:', mycell.path
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_suppyrrs', 'Vm', sim.data)

        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=5e-3, firstWidth=100e-3)

        sim.schedule()
        if mycell.has_cycle():
            print "WARNING!! CYCLE PRESENT IN CICRUIT."
        t1 = datetime.now()
        sim.run(10e-3)
        t2 = datetime.now()
        delta = t2 - t1
        print 'simulation time: ', delta.seconds + 1e-6 * delta.microseconds
        sim.dump_data('data')
        mycell.dump_cell('suppyrrs.txt')
        print 'soma:', 'Ra =', mycell.soma.Ra, 'Rm =', mycell.soma.Rm, 'Cm =', mycell.soma.Cm, 'Em =', mycell.soma.Em, 'initVm =', mycell.soma.initVm
        print 'dend:', 'Ra =', mycell.comp[2].Ra, 'Rm =', mycell.comp[2].Rm, 'Cm =', mycell.comp[2].Cm, 'Em =', mycell.comp[2].Em, 'initVm =', mycell.comp[2].initVm
        mus_vm = pylab.array(vm_table) * 1e3
        pylab.plot(mus_vm, 'r-', label='mus')
        pylab.legend()
        pylab.show()
        
        
# test main --
from simulation import Simulation
import pylab

if __name__ == "__main__":
    SupPyrRS.test_single_cell()
    

# 
# suppyrrs.py ends here
