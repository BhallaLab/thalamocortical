# test_singlecomp.py --- 
# 
# Filename: test_singlecomp.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Jul 17 17:50:19 2012 (+0530)
# Version: 
# Last-Updated: Thu Jul 19 17:00:15 2012 (+0530)
#           By: subha
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

import sys
sys.path.append('/data/subha/chamcham_moose/python')

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
    # 'NaPF_SS':  1.5,
    # 'KDR_FS':   1000.0,
    # 'KC_FAST':  100.0,
    # 'KA':       300.0,
    # 'KM':       37.5,
    # 'K2':       1.0,
    # 'KAHP_SLOWER':      1.0,
    # 'CaL':      5.0,
    # 'CaT_A':    1.0,
    # 'AR':       2.5
}
# channels = {'NaF2': 'NaF2_SS', 'NaPF_SS': 'NaPF_SS', 'KDR_FS': 'KDR_FS', \
#                 'KA': 'KA', 'K2': 'K2', 'KM': 'KM', 'KC_FAST': 'KC_FAST', \
#                 'KAHP_SLOWER': 'KAHP_SLOWER', \
#                 'CaL': 'CaL', 'CaT_A': 'CaT_A', 'AR': 'AR'}

ENa = 50e-3
EK = -100e-3
ECa = 125e-3
Em = -65e-3
EAR = -40e-3

simtime = 350e-3
simdt = 0.25e-4
plotdt = 0.25e-4

channels_inited = False
def init_channels():
    global channels_inited
    global channels
    if channels_inited:
        return
        
    lib = moose.Neutral('/library')
    channel_lib = {}
    channel = None
    for channel_class in channel_density.keys():
        if config.context.exists('/library/' + channel_class):
            channel = moose.HHChannel(channel_class, lib)
        else:
            class_obj = eval(channel_class)
            if channel_class == 'NaF2':
                channel = class_obj(channel_class, lib, shift=-2.5e-3)
            else:
                channel = class_obj(channel_class, lib)
            channel.X = 0.0
        channel_lib[channel_class] = channel
    channels_inited = True
    return channel_lib
            
if __name__ == '__main__':
    sim = Simulation('test_ss_soma')
    soma = MyCompartment('soma', sim.model)
    soma.length = 20e-6
    soma.diameter = 2e-6 * 7.5
    soma.setSpecificCm(9e-3)
    soma.setSpecificRm(5.0)
    soma.setSpecificRa(1.0)
    soma.Em = -65e-3
    soma.initVm = -65e-3
    channel_lib = init_channels()
    gk = {}
    for channel, density in channel_density.items():
        chan = channel_lib[channel]
        new_chan = moose.HHChannel(chan, chan.name, soma)
        chan = soma.insertChannel(new_chan, density)
        chan.X = 0.0
        gk[channel] = moose.Table(channel, sim.data)
        gk[channel].stepMode = 3
        gk[channel].connect('inputRequest', chan, 'Gk')
        print chan.name, chan.Gbar
	if channel.startswith('K'):
	    chan.Ek = EK
	elif channel.startswith('Na'):
	    chan.Ek = ENa
	elif channel.startswith('Ca'):
	    chan.Ek = ECa
	elif channel.startswith('AR'):
	    chan.Ek = EAR
#   	    chan.X = 0.25
	else:
	    print 'Error: unknown channel', channel
    cad = moose.CaConc('CaPool', soma)
    cad.B = (5.2e-6 / 2e-10) / soma.sarea()
    cad.tau = 50e-3
    vm_table = soma.insertRecorder('Vm', 'Vm', sim.data)
    soma.insertPulseGen('pulsegen', sim.model, firstLevel=3e-10, firstDelay=1e9, firstWidth=1e3)
    print 'Rm:', soma.Rm
    print 'Cm:', soma.Cm,
    print 'Em:', soma.Em
    print 'initVm:', soma.initVm
    gk_naf2_table = moose.Table('Gk_NaF2', sim.data)
    gk_naf2_table.stepMode = 3
    chan = moose.HHChannel(soma.path + '/NaF2')
    print chan.Gbar, chan.Ek
    chan.connect('Gk', gk_naf2_table, 'inputRequest')
    sim.schedule(simdt=simdt, plotdt=plotdt)
    print 'Soma fields'
    for field in dir(soma):
        x = getattr(soma, field)
        if isinstance(x, float):
            print field, x
    chan = moose.HHChannel(soma.path + '/NaF2')
    print 'NaF2 fields'
    for field in dir(chan):
        x = getattr(chan, field)
        if isinstance(x, float):
            print field, x
    # for channel in soma.channels:
    #     channel.useClock(0)
    sim.run(simtime)
    for i in range(len(gk_naf2_table)):
        gk_naf2_table[i] = gk_naf2_table[i] / soma.sarea()
    sim.dump_data('data')
    pylab.subplot(111)
    tseries = pylab.linspace(0, simtime, len(vm_table))
    pylab.plot(tseries*1e3, pylab.array(vm_table) * 1e3, 'x', label='oldmus')
    # pylab.subplot(312)
    newmus = pylab.loadtxt('/home/subha/src/dh_branch/Demos/traub_2005/py/data/singlecomp_Vm.dat')
    pylab.plot(newmus[:,0], newmus[:,1]*1e3, 'x', label='newmus')
    nrn = pylab.loadtxt('/home/subha/src/dh_branch/Demos/traub_2005/nrn/data/singlecomp_Vm.plot')
    pylab.plot(nrn[:,0], nrn[:,1], 'x', label='nrn')
    # for key, value in gk.items():
    #     pylab.plot(value, label=key)
    pylab.legend()
    pylab.show()



# 
# test_singlecomp.py ends here
