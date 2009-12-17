# spinystellate.py --- 
# 
# Filename: spinystellate.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Tue Sep 29 11:43:22 2009 (+0530)
# Version: 
# Last-Updated: Thu Nov 12 14:38:37 2009 (+0530)
#           By: subhasis ray
#     Update #: 149
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
import moose
import config
from cell import *
from capool import CaPool
# from cellview import MyCellView

class SpinyStellate(TraubCell):
    ENa = 50e-3
    EK = -100e-3
    EAR = -40e-3
    ECa = 100e-3
    prototype = TraubCell.read_proto("SpinyStellate.p", "SpinyStellate")
    def __init__(self, *args):
	TraubCell.__init__(self, *args)
        self.ar = moose.HHChannel('/model/SpinyStellate/comp_1/AR')

    def _topology(self):
	self.presyn = 57
        self.level[1].add(self.comp[1])
        for i in range(2, 42, 13):
            self.level[2].add(self.comp[i])
        for i in range(3, 43, 13):
            self.level[3].add(self.comp[i])
            self.level[3].add(self.comp[i+1])
        for i in range(5, 45, 13):
            self.level[4].add(self.comp[i])
            self.level[4].add(self.comp[i+1])
            self.level[4].add(self.comp[i+2])
        for i in range(8, 48, 13):
            self.level[5].add(self.comp[i])
            self.level[5].add(self.comp[i+1])
            self.level[5].add(self.comp[i+2])
        for i in range(11, 51, 13):
            self.level[6].add(self.comp[i])
            self.level[7].add(self.comp[i+1])
            self.level[8].add(self.comp[i+2])
            self.level[9].add(self.comp[i+3])

        for i in range(54, 60):
            self.level[0].add(self.comp[i])
            
        # Skipping the categorizatioon into levels for the time being

    def _setup_passive(self):
	for comp in self.comp[1:]:
	    comp.initVm = -65e-3

    def _setup_channels(self):
        """Set up connection between CaPool, Ca channels, Ca dependnet channels."""
        for comp in self.comp[1:]:
            ca_pool = None
            ca_dep_chans = []
            ca_chans = []
            for child in comp.children():
                obj = moose.Neutral(child)
                if obj.name == 'CaPool':
                    ca_pool = moose.CaConc(child)
                    ca_pool.tau = 20e-3
                elif obj.className == 'HHChannel':
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
        sim = Simulation()
        mycell = SpinyStellate(SpinyStellate.prototype, sim.model.path + "/SpinyStellate")
        mycell.soma.x0 = 0.0
        mycell.soma.y0 = 0.0
        mycell.soma.z0 = 0.0
        mycell.soma.x = 0.0
        mycell.soma.y = 0.0
        mycell.soma.z = mycell.soma.length
        mycellview = MyCellView(mycell)
        print 'Created cell:', mycell.path
        for neighbour in mycell.soma.neighbours('raxial'):
            print 'RAXIAL', neighbours.path()
        for neighbour in mycell.soma.neighbours('axial'):
            print 'AXIAL', neighbour.path()
        vm_table = mycell.comp[mycell.presyn].insertRecorder('Vm_spinstell', 'Vm', sim.data)
        pulsegen = mycell.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=50e-3, firstWidth=50e-3)
#         pulsegen1 = mycell.soma.insertPulseGen('pulsegen1', sim.model, firstLevel=3e-7, firstDelay=150e-3, firstWidth=10e-3)

        sim.schedule()
        if mycell.has_cycle():
            print "WARNING!! CYCLE PRESENT IN CICRUIT."
        t1 = datetime.now()
        sim.run(200e-3)
        t2 = datetime.now()
        delta = t2 - t1
        print 'simulation time: ', delta.seconds + 1e-6 * delta.microseconds
        sim.dump_data('data')
        mycell.dump_cell('spinstell.txt')
        
        mus_vm = pylab.array(vm_table) * 1e3
        nrn_vm = pylab.loadtxt('../nrn/mydata/Vm_spinstell.plot')
        nrn_t = nrn_vm[:, 0]
        mus_t = linspace(0, nrn_t[-1], len(mus_vm))
        nrn_vm = nrn_vm[:, 1]
        nrn_ca = pylab.loadtxt('../nrn/mydata/Ca_spinstell.plot')
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
#     call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_spinstell.hoc'], cwd='../nrn')
    SpinyStellate.test_single_cell()
    


# 
# spinystellate.py ends here
