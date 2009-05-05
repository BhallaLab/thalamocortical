# test_ss_soma.py --- 
# 
# Filename: test_ss_soma.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sun May  3 12:52:23 2009 (+0530)
# Version: 
# Last-Updated: Tue May  5 17:11:33 2009 (+0530)
#           By: subhasis ray
#     Update #: 129
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
import pylab

from simulation import Simulation
from kchans import *
from nachans import *
from cachans import *
from capool import *
from archan import *

from compartment import MyCompartment


channel_density = {
    'NaF2':     1500.0,
    'NaPF_SS':  1.5,
    'KDR_FS':   1000.0,
    'KC_FAST':  100.0,
    'KA':       300.0,
    'KM':       37.5,
    'K2':       1.0,
    'KAHP_SLOWER':      1.0,
    'CaL':      5.0,
    'CaT_A':    1.0,
    'AR':       2.5
}
channels = {'NaF2': 'NaF2_SS', 'NaPF_SS': 'NaPF_SS', 'KDR_FS': 'KDR_FS', \
                'KA': 'KA', 'K2': 'K2', 'KM': 'KM', 'KC_FAST': 'KC_FAST', \
                'KAHP_SLOWER': 'KAHP_SLOWER', \
                'CaL': 'CaL', 'CaT_A': 'CaT_A', 'AR': 'AR'}

ENa = 50e-3
EK = -100e-3
ECa = 125e-3
Em = -65e-3
EAR = -40e-3
channels_inited = False
def init_channels():
    global channels_inited
    global channels
    if channels_inited:
        return
        
    lib = moose.Neutral('/library')
    channel_lib = {}
    for channel_class, channel_name in channels.items():
        if config.context.exists('/library/' + channel_name):
            channel = moose.HHChannel(channel_name, lib)
        else:
            class_obj = eval(channel_class)
            if channel_class == 'NaF2':
                channel = class_obj(channel_name, lib, shift=-2.5e-3)
            else:
                channel = class_obj(channel_name, lib)
        channel_lib[channel_class] = channel
    channels_inited = True
    return channel_lib
            
if __name__ == '__main__':
    sim = Simulation()
    soma = MyCompartment('soma', sim.model)
    soma.length = 40e-6
    soma.diameter = 2e-6 * 7.5
    soma.setSpecificCm(9e-3)
    soma.setSpecificRm(5.0)
    soma.setSpecificRa(1.0)
    soma.Em = -65e-3
    soma.initVm = -65e-3
    channel_lib = init_channels()
    for channel, density in channel_density.items():
        chan = channel_lib[channel]
	chan = soma.insertChannel(chan, density)
	if channel.startswith('K'):
	    chan.Ek = EK
	elif channel.startswith('Na'):
	    chan.Ek = ENa
	elif channel.startswith('Ca'):
	    chan.Ek = ECa
	elif channel.startswith('AR'):
	    chan.Ek = EAR
  	    chan.X = 0.25
	else:
	    print 'Error: unknown channel', channel
	    
    soma.insertCaPool(5.2e-6 / 2e-10, 50e-3)
    vm_table = soma.insertRecorder('Vm', sim.data)
    soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=20.0e-3, firstWidth=1e3)
    sim.schedule()
    config.context.useClock(0, sim.model.path + '/##')
    soma.useClock(1, 'init')
    for channel in soma.channels:
        channel.useClock(0)
    sim.run(50e-3)
    sim.dump_data('data')
    pylab.plot(vm_table)
    pylab.show()



# 
# test_ss_soma.py ends here
