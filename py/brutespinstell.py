# brutespinstell.py --- 
# 
# Filename: brutespinstell.py
# Description: This is an unelegant version of spiny stellate cells
# Author: subhasis ray
# Maintainer: 
# Created: Fri May  8 11:24:30 2009 (+0530)
# Version: 
# Last-Updated: Fri May  8 21:26:42 2009 (+0530)
#           By: subhasis ray
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

from collections import deque, defaultdict
from datetime import datetime
import moose

from kchans import *
from nachans import *
from cachans import *
from capool import *
from archan import *

from compartment import MyCompartment

class SpinyStellate(moose.Neutral):
    ENa = 50e-3
    EK = -100e-3
    ECa = 125e-3
    Em = -65e-3
    EAR = -40e-3
    conductance ={
	0: {
	    'NaF2':   0.4,
	    'KDR_FS':   0.4,
	    'KA':   0.002,
	    'K2':   0.0001
	    },
	1: {
	    'NaF2':   0.15,
	    'NaPF_SS':   0.00015,
	    'KDR_FS':   0.1,
	    'KC_FAST':   0.01,
	    'KA':   0.03,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	2: {
	    'NaF2':   0.075,
	    'NaPF_SS':   7.5E-05,
	    'KDR_FS':   0.075,
	    'KC_FAST':   0.01,
	    'KA':   0.03,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	3: {
	    'NaF2':   0.075,
	    'NaPF_SS':   7.5E-05,
	    'KDR_FS':   0.075,
	    'KC_FAST':   0.01,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	4: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KC_FAST':   0.01,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	5: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	6: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.0005,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	7: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.003,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	8: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.003,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    },
	9: {
	    'NaF2':   0.005,
	    'NaPF_SS':   5.E-06,
	    'KA':   0.002,
	    'KM':   0.00375,
	    'K2':   0.0001,
	    'KAHP_SLOWER':   0.0001,
	    'CaL':   0.003,
	    'CaT_A':   0.0001,
	    'AR':   0.00025
	    }
	}

    channels = {'NaF2': 'NaF2_SS', 
                'NaPF_SS': 'NaPF_SS',
                'KDR_FS': 'KDR_FS',
                'KA': 'KA', 
                'K2': 'K2', 
                'KM': 'KM', 
                'KC_FAST': 'KC_FAST', 
                'KAHP_SLOWER': 'KAHP_SLOWER',
                'CaL': 'CaL', 
                'CaT_A': 'CaT_A', 
                'AR': 'AR'}

    def __init__(self, *args):
	moose.Neutral.__init__(self, *args)
        self.channel_lib = {}
        for channel_class, channel_name in SpinyStellate.channels.items():
            channel = None
            if config.context.exists('/library/' + channel_name):
                channel = moose.HHChannel(channel_name, config.lib)
            else:
                class_obj = eval(channel_class)
                if channel_class == 'NaF2':
                    channel = class_obj(channel_name, config.lib, shift=0.0)
                else:
                    channel = class_obj(channel_name, config.lib)
            self.channel_lib[channel_class] = channel

	comp = []
	dendrites = set()
	level = defaultdict(set)
	axon = []
	self.presyn = 57
	for ii in range(60):
	    comp.append(MyCompartment('comp_' + str(ii), self))
	comp[1].connect('raxial', comp[ 54], 'axial')
	comp[1].connect('raxial', comp[ 2], 'axial') 
	comp[1].connect('raxial', comp[ 15], 'axial')
	comp[1].connect('raxial', comp[ 28], 'axial')
	comp[1].connect('raxial', comp[ 41], 'axial')
	comp[2].connect('raxial', comp[ 3], 'axial') 
	comp[2].connect('raxial', comp[ 4], 'axial') 
	comp[3].connect('raxial', comp[ 4], 'axial') 
	comp[3].connect('raxial', comp[ 5], 'axial') 
	comp[3].connect('raxial', comp[ 6], 'axial') 
	comp[4].connect('raxial', comp[ 7], 'axial') 
	comp[5].connect('raxial', comp[ 6], 'axial') 
	comp[5].connect('raxial', comp[ 8], 'axial') 
	comp[6].connect('raxial', comp[ 9], 'axial') 
	comp[7].connect('raxial', comp[ 10], 'axial')
	comp[8].connect('raxial', comp[ 11], 'axial')
	comp[11].connect('raxial', comp[  12], 'axial')
	comp[12].connect('raxial', comp[  13], 'axial')
	comp[13].connect('raxial', comp[  14], 'axial')
	comp[15].connect('raxial', comp[  16], 'axial')
	comp[15].connect('raxial', comp[  17], 'axial')
	comp[16].connect('raxial', comp[  17], 'axial')
	comp[16].connect('raxial', comp[  18], 'axial')
	comp[16].connect('raxial', comp[  19], 'axial')
	comp[17].connect('raxial', comp[  20], 'axial')
	comp[18].connect('raxial', comp[  19], 'axial')
	comp[18].connect('raxial', comp[  21], 'axial')
	comp[19].connect('raxial', comp[  22], 'axial')
	comp[20].connect('raxial', comp[  23], 'axial')
	comp[21].connect('raxial', comp[  24], 'axial')
	comp[24].connect('raxial', comp[  25], 'axial')
	comp[25].connect('raxial', comp[  26], 'axial')
	comp[26].connect('raxial', comp[  27], 'axial')
	comp[28].connect('raxial', comp[  29], 'axial')
	comp[28].connect('raxial', comp[  30], 'axial')
	comp[29].connect('raxial', comp[  30], 'axial')
	comp[29].connect('raxial', comp[  31], 'axial')
	comp[29].connect('raxial', comp[  32], 'axial')
	comp[30].connect('raxial', comp[  33], 'axial')
	comp[31].connect('raxial', comp[  32], 'axial')
	comp[31].connect('raxial', comp[  34], 'axial')
	comp[32].connect('raxial', comp[  35], 'axial')
	comp[33].connect('raxial', comp[  36], 'axial')
	comp[34].connect('raxial', comp[  37], 'axial')
	comp[37].connect('raxial', comp[  38], 'axial')
	comp[38].connect('raxial', comp[  39], 'axial')
	comp[39].connect('raxial', comp[  40], 'axial')
	comp[41].connect('raxial', comp[  42], 'axial')
	comp[41].connect('raxial', comp[  43], 'axial')
	comp[42].connect('raxial', comp[  43], 'axial')
	comp[42].connect('raxial', comp[  44], 'axial')
	comp[42].connect('raxial', comp[  45], 'axial')
	comp[43].connect('raxial', comp[  46], 'axial')
	comp[44].connect('raxial', comp[  45], 'axial')
	comp[44].connect('raxial', comp[  47], 'axial')
	comp[45].connect('raxial', comp[  48], 'axial')
	comp[46].connect('raxial', comp[  49], 'axial')
	comp[47].connect('raxial', comp[  50], 'axial')
	comp[50].connect('raxial', comp[  51], 'axial')
	comp[51].connect('raxial', comp[  52], 'axial')
	comp[52].connect('raxial', comp[  53], 'axial')
	comp[54].connect('raxial', comp[  55], 'axial')
	comp[55].connect('raxial', comp[  56], 'axial')
	comp[55].connect('raxial', comp[  58], 'axial')
	comp[56].connect('raxial', comp[  57], 'axial')
	comp[56].connect('raxial', comp[  58], 'axial')
	comp[58].connect('raxial', comp[  59], 'axial')

	comp[ 1].diameter = 2 * 7.5 
	comp[ 2].diameter = 2 * 1.06 
	comp[ 3].diameter = 2 * 0.666666667 
	comp[ 4].diameter = 2 * 0.666666667 
	comp[ 5].diameter = 2 * 0.418972332 
	comp[ 6].diameter = 2 * 0.418972332 
	comp[ 7].diameter = 2 * 0.666666667 
	comp[ 8].diameter = 2 * 0.418972332 
	comp[ 9].diameter = 2 * 0.418972332 
	comp[ 10].diameter = 2 * 0.666666667 
	comp[ 11].diameter = 2 * 0.418972332 
	comp[ 12].diameter = 2 * 0.418972332 
	comp[ 13].diameter = 2 * 0.418972332 
	comp[ 14].diameter = 2 * 0.418972332 
	comp[ 15].diameter = 2 * 1.06 
	comp[ 16].diameter = 2 * 0.666666667 
	comp[ 17].diameter = 2 * 0.666666667 
	comp[ 18].diameter = 2 * 0.418972332 
	comp[ 19].diameter = 2 * 0.418972332 
	comp[ 20].diameter = 2 * 0.666666667 
	comp[ 21].diameter = 2 * 0.418972332 
	comp[ 22].diameter = 2 * 0.418972332 
	comp[ 23].diameter = 2 * 0.666666667 
	comp[ 24].diameter = 2 * 0.418972332 
	comp[ 25].diameter = 2 * 0.418972332 
	comp[ 26].diameter = 2 * 0.418972332 
	comp[ 27].diameter = 2 * 0.418972332 
	comp[ 28].diameter = 2 * 1.06 
	comp[ 29].diameter = 2 * 0.666666667 
	comp[ 30].diameter = 2 * 0.666666667 
	comp[ 31].diameter = 2 * 0.418972332 
	comp[ 32].diameter = 2 * 0.418972332 
	comp[ 33].diameter = 2 * 0.666666667 
	comp[ 34].diameter = 2 * 0.418972332 
	comp[ 35].diameter = 2 * 0.418972332 
	comp[ 36].diameter = 2 * 0.666666667 
	comp[ 37].diameter = 2 * 0.418972332 
	comp[ 38].diameter = 2 * 0.418972332 
	comp[ 39].diameter = 2 * 0.418972332 
	comp[ 40].diameter = 2 * 0.418972332 
	comp[ 41].diameter = 2 * 1.06 
	comp[ 42].diameter = 2 * 0.666666667 
	comp[ 43].diameter = 2 * 0.666666667 
	comp[ 44].diameter = 2 * 0.418972332 
	comp[ 45].diameter = 2 * 0.418972332 
	comp[ 46].diameter = 2 * 0.666666667 
	comp[ 47].diameter = 2 * 0.418972332 
	comp[ 48].diameter = 2 * 0.418972332 
	comp[ 49].diameter = 2 * 0.666666667 
	comp[ 50].diameter = 2 * 0.418972332 
	comp[ 51].diameter = 2 * 0.418972332 
	comp[ 52].diameter = 2 * 0.418972332 
	comp[ 53].diameter = 2 * 0.418972332 
	comp[ 54].diameter = 2 * 0.7 
	comp[ 55].diameter = 2 * 0.6 
	comp[ 56].diameter = 2 * 0.5 
	comp[ 57].diameter = 2 * 0.5 
	comp[ 58].diameter = 2 * 0.5 
	comp[ 59].diameter = 2 * 0.5 
	comp[ 1].length = 20. 
	comp[ 2].length = 40. 
	comp[ 3].length = 40. 
	comp[ 4].length = 40. 
	comp[ 5].length = 40. 
	comp[ 6].length = 40. 
	comp[ 7].length = 40. 
	comp[ 8].length = 40. 
	comp[ 9].length = 40. 
	comp[ 10].length = 40.
	comp[ 11].length = 40.
	comp[ 12].length = 40.
	comp[ 13].length = 40.
	comp[ 14].length = 40.
	comp[ 15].length = 40.
	comp[ 16].length = 40.
	comp[ 17].length = 40.
	comp[ 18].length = 40.
	comp[ 19].length = 40.
	comp[ 20].length = 40.
	comp[ 21].length = 40.
	comp[ 22].length = 40.
	comp[ 23].length = 40.
	comp[ 24].length = 40.
	comp[ 25].length = 40.
	comp[ 26].length = 40.
	comp[ 27].length = 40.
	comp[ 28].length = 40.
	comp[ 29].length = 40.
	comp[ 30].length = 40.
	comp[ 31].length = 40.
	comp[ 32].length = 40.
	comp[ 33].length = 40.
	comp[ 34].length = 40.
	comp[ 35].length = 40.
	comp[ 36].length = 40.
	comp[ 37].length = 40.
	comp[ 38].length = 40.
	comp[ 39].length = 40.
	comp[ 40].length = 40.
	comp[ 41].length = 40.
	comp[ 42].length = 40.
	comp[ 43].length = 40.
	comp[ 44].length = 40.
	comp[ 45].length = 40.
	comp[ 46].length = 40.
	comp[ 47].length = 40.
	comp[ 48].length = 40.
	comp[ 49].length = 40.
	comp[ 50].length = 40.
	comp[ 51].length = 40.
	comp[ 52].length = 40.
	comp[ 53].length = 40.
	comp[ 54].length = 50.
	comp[ 55].length = 50.
	comp[ 56].length = 50.
	comp[ 57].length = 50.
	comp[ 58].length = 50.
	comp[ 59].length = 50.
	
	level[ 1].add( comp[ 1]) 
	level[ 2].add( comp[ 2]) 
	level[ 3].add( comp[ 3]) 
	level[ 3].add( comp[ 4]) 
	level[ 4].add( comp[ 5]) 
	level[ 4].add( comp[ 6]) 
	level[ 4].add( comp[ 7]) 
	level[ 5].add( comp[ 8]) 
	level[ 5].add( comp[ 9]) 
	level[ 5].add( comp[ 10])
	level[ 6].add( comp[ 11])
	level[ 7].add( comp[ 12])
	level[ 8].add( comp[ 13])
	level[ 9].add( comp[ 14])
	level[ 2].add( comp[ 15])
	level[ 3].add( comp[ 16])
	level[ 3].add( comp[ 17])
	level[ 4].add( comp[ 18])
	level[ 4].add( comp[ 19])
	level[ 4].add( comp[ 20])
	level[ 5].add( comp[ 21])
	level[ 5].add( comp[ 22])
	level[ 5].add( comp[ 23])
	level[ 6].add( comp[ 24])
	level[ 7].add( comp[ 25])
	level[ 8].add( comp[ 26])
	level[ 9].add( comp[ 27])
	level[ 2].add( comp[ 28])
	level[ 3].add( comp[ 29])
	level[ 3].add( comp[ 30])
	level[ 4].add( comp[ 31])
	level[ 4].add( comp[ 32])
	level[ 4].add( comp[ 33])
	level[ 5].add( comp[ 34])
	level[ 5].add( comp[ 35])
	level[ 5].add( comp[ 36])
	level[ 6].add( comp[ 37])
	level[ 7].add( comp[ 38])
	level[ 8].add( comp[ 39])
	level[ 9].add( comp[ 40])
	level[ 2].add( comp[ 41])
	level[ 3].add( comp[ 42])
	level[ 3].add( comp[ 43])
	level[ 4].add( comp[ 44])
	level[ 4].add( comp[ 45])
	level[ 4].add( comp[ 46])
	level[ 5].add( comp[ 47])
	level[ 5].add( comp[ 48])
	level[ 5].add( comp[ 49])
	level[ 6].add( comp[ 50])
	level[ 7].add( comp[ 51])
	level[ 8].add( comp[ 52])
	level[ 9].add( comp[ 53])
	level[ 0].add( comp[ 54])
	level[ 0].add( comp[ 55])
	level[ 0].add( comp[ 56])
	level[ 0].add( comp[ 57])
	level[ 0].add( comp[ 58])
	level[ 0].add( comp[ 59])
	
	for ii in range(2, len(level)):
	    dendrites |= level[ii]
	self.level = level
	self.comp = comp
	self.dendrites = dendrites
	self.soma = comp[1]
	for compartment in comp[1:]:
	    compartment.length *= 1e-6
	    compartment.diameter *= 1e-6
	    compartment.setSpecificCm(9e-3)


        t1 = datetime.now()
	    
	for ii in range(0, len(level)):
	    comp_set = level[ii]
	    conductances = SpinyStellate.conductance[ii]
	    mult = 1e4
	    if ii > 1:
		mult *= 2.0
	    for comp in comp_set:
		for channel_name, density in conductances.items():
		    channel = moose.HHChannel(self.channel_lib[channel_name], 
					  channel_name, comp)
		    comp.insertChannel(channel, specificGbar=mult *  density)
                    if channel_name.startswith('K'):
                        channel.Ek = SpinyStellate.EK
                    elif channel_name.startswith('Na'):
                        channel.X = 0.0
                        channel.Ek = SpinyStellate.ENa
                    elif channel_name.startswith('Ca'):
                        channel.Ek = SpinyStellate.ECa
                    elif channel_name.startswith('AR'):
                        channel.Ek = SpinyStellate.EAR
                        channel.X = 0.0
                    else:
                        print 'ERROR: Unknown channel type:', channel
	for compartment in self.dendrites:
	    print compartment.name, compartment.length, compartment.diameter
	    compartment.setSpecificRm(5.0/2)
	    compartment.setSpecificRa(2.5)
	    compartment.Cm *= 2.0
	    compartment.insertCaPool(5.2e-6 / 2e-10, 20e-3)

	self.soma.setSpecificRm(5.0)
	self.soma.setSpecificRa(2.5)
        self.soma.insertCaPool(5.2e-6 / 2e-10, 50e-3)

	for compartment in self.level[0]: # axonal comps
	    compartment.setSpecificRm(0.1)
	    compartment.setSpecificRa(1.0)

        t2 = datetime.now()
        delta = t2 - t1
        print 'insert channels: ', delta.seconds + 1e-6 * delta.microseconds
		    

import pylab
from simulation import Simulation
if __name__ == '__main__':
    sim = Simulation()
    s = SpinyStellate('cell', sim.model)
    vm_table = s.soma.insertRecorder('Vm', sim.data)
    pulsegen = s.soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=0.0, firstWidth=50e-3)
    sim.schedule()
    t1 = datetime.now()
    sim.run(50e-3)
    t2 = datetime.now()
    delta = t2 - t1
    print 'simulation time: ', delta.seconds + 1e-6 * delta.microseconds
    sim.dump_data('data')
    pylab.plot(vm_table)
    pylab.show()

# 
# brutespinstell.py ends here
