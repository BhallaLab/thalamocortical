# testutils.py --- 
# 
# Filename: testutils.py
# Description: 
# Author: 
# Maintainer: 
# Created: Sat May 26 12:02:18 2012 (+0530)
# Version: 
# Last-Updated: Thu May 31 15:40:58 2012 (+0530)
#           By: subha
#     Update #: 35
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

import moose

def make_testcomp(containerpath):
    comp = moose.Compartment('%s/testcomp' % (containerpath))
    comp.Em = -65e-3
    comp.initVm = -65e-3
    comp.Cm = 1e-12
    comp.Rm = 1e9
    comp.Ra = 1e5
    return comp

def make_pulsegen(containerpath):
    pulsegen = moose.PulseGen('%s/testpulse' % (containerpath))
    pulsegen.firstLevel = 1e-12
    pulsegen.firstDelay = 50e-3
    pulsegen.firstWidth = 100e-3
    pulsegen.secondLevel = -1e-12
    pulsegen.secondDelay = 150e-3
    pulsegen.secondWidth = 100e-3
    return pulsegen

def setup_single_compartment(container_path, channel_proto, Gbar=1e-9):
    comp = make_testcomp(container_path)
    moose.context.copy(channel_proto.id, comp.id, channel_proto.name)
    channel = moose.HHChannel('%s/%s' %(comp.path, channel_proto.name))
    channel.connect('channel', comp, 'channel')
    channel.Gbar = Gbar
    pulsegen = make_pulsegen(container_path)
    pulsegen.connect('outputSrc', comp, 'injectMsg')
    vm_table = moose.Table('%s/Vm' % (container_path))
    vm_table.connect('inputRequest', comp, 'Vm')
    vm_table.stepMode = 3
    gk_table = moose.Table('%s/Gk' % (container_path))
    gk_table.connect('inputRequest', channel, 'Gk')
    gk_table.stepMode = 3
    moose.context.setClock(0, 1e-5)
    moose.context.setClock(1, 1e-5)
    moose.context.useClock(0, '%s/##[TYPE=Compartment]' % (container_path), 'init')
    moose.context.useClock(1, '%s/##' % (container_path), 'process')
    return {'compartment': comp,
            'stimulus': pulsegen,
            'channel': channel,
            'Vm': vm_table,
            'Gk': gk_table}



# 
# testutils.py ends here
