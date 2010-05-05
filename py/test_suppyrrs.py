# test_suppyrrs.py --- 
# 
# Filename: test_suppyrrs.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Sep 18 16:38:41 2009 (+0530)
# Version: 
# Last-Updated: Sat Sep 19 17:13:21 2009 (+0530)
#           By: subhasis ray
#     Update #: 67
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# Test for single compartment model ... just NaF and CaL interaction
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
call(['/home/subha/neuron/nrn/x86_64/bin/nrngui', 'test_suppyrRS_1comp.hoc'], cwd='../nrn')
from datetime import datetime
import moose
import config
from cell import *
from capool import CaPool
from compartment import MyCompartment
from simulation import Simulation
sim = Simulation()
cell = moose.Cell('cell', sim.model)
soma = MyCompartment('soma', cell)
soma.diameter = 16e-6
soma.length = 15e-6
soma.Em = -70e-3
soma.initVm = -65e-3
soma.setSpecificRm(5.0)
soma.setSpecificRa(2.5)
soma.setSpecificCm(9e-3)
soma.insertChannel('NaF', specificGbar=1875.0, Ek=50e-3, shift=-3.5e-3)
soma.insertChannel('CaL', specificGbar=10.0, Ek=125e-3)
soma.insertCaPool(2600000.0, 100e-3)
vmTable = soma.insertRecorder('Vm', 'Vm', sim.data)
caTable = moose.Table('ca', sim.data)
caTable.stepMode = 3
soma.ca_pool.connect('Ca', caTable, 'inputRequest')
sim.schedule()
sim.run(200e-3)
sim.dump_data('data')

from pylab import *
mus_Ca = array(caTable)
mus_data = loadtxt('data/Vm.plot')
nrn_data = loadtxt('../nrn/mydata/Vm.plot')
nrn_Vm = nrn_data[:,1]
nrn_t = nrn_data[:,0]
nrn_data = loadtxt('../nrn/mydata/Ca.plot')
nrn_Ca = nrn_data[:,1]
plot(nrn_t, nrn_Ca, 'b-', label='nrn')
mus_Vm = mus_data
mus_t = linspace(0, nrn_t[-1], len(mus_Vm))
plot(mus_t, mus_Ca, 'y-',label='mus')
legend()
show()



# 
# test_suppyrrs.py ends here
