# test_cat_a.py --- 
# 
# Filename: test_cat_a.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Thu Sep 24 18:23:34 2009 (+0530)
# Version: 
# Last-Updated: Fri Sep 25 17:07:24 2009 (+0530)
#           By: subhasis ray
#     Update #: 48
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

from subprocess import call
from pylab import *
from moose import *
from cachans import *
from simulation import *
from compartment import *
from config import context
if __name__ == '__main__':
    sim = Simulation()
    comp = MyCompartment('comp')
    comp.length = 20e-6
    comp.diameter = 15e-6
    comp.setSpecificRm(2.5)
    comp.setSpecificRa(2.0)
    comp.setSpecificCm(0.01)
    chan = comp.insertChannel('CaT_A', specificGbar=5)
    comp.Em = -65e-3
    comp.initVm = -65e-3
    comp.insertPulseGen('pulsegen', sim.model, firstLevel=.1e-9 , firstWidth=20e-3, firstDelay=20e-3)
    comp.insertRecorder('Vm', 'Vm', sim.data)
    sim.schedule()
    context.reset()
    sim.run(50e-3)
    sim.dump_data('data')
#     call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_cat_a.hoc'], cwd='../nrn')
    mus = loadtxt('data/Vm.plot')
    nrn = loadtxt('../nrn/mydata/Vm.plot')
    plot(linspace(0.0, 50.0, len(mus)), mus*1e3, 'r-', label='mus')
    plot(nrn[:,0], nrn[:,1], 'b-.', label='nrn')
    legend()
    show()
    

# 
# test_cat_a.py ends here
