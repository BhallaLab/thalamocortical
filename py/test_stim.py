# test_stim.py --- 
# 
# Filename: test_stim.py
# Description: 
# Author: 
# Maintainer: 
# Created: Thu Jan  5 16:59:17 2012 (+0530)
# Version: 
# Last-Updated: Thu Jan  5 18:29:41 2012 (+0530)
#           By: subha
#     Update #: 28
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
import pylab
import numpy
def test_stim():
    stim_gate = moose.PulseGen('gate')
    stim_gate.trigMode = moose.FREE_RUN
    stim_gate.firstLevel = 1.0
    stim_gate.firstDelay = 1.0
    stim_gate.firstWidth = 1e9
    table_gate = moose.Table('tab_gate')
    table_gate.stepMode = 3
    table_gate.connect('inputRequest', stim_gate, 'output')
    stim_bg = moose.PulseGen('bg')
    stim_bg.trigMode = moose.EXT_GATE
    stim_bg.firstDelay = 0.5
    stim_bg.firstLevel = 1.0
    stim_bg.firstWidth = 60e-6
    stim_bg.secondLevel = 1.0
    stim_bg.secondWidth = 60e-6
    stim_bg.secondDelay = 10e-3
    tab_bg = moose.Table('tab_bg')
    tab_bg.stepMode = 3
    tab_bg.connect('inputRequest', stim_bg, 'output')
    stim_gate.connect('outputSrc', stim_bg, 'input')
    # stim_probe = moose.PulseGen('probe')
    # stim_probe.trigMode = moose.EXT_GATE
    # stim_probe.firstLevel = 0.5
    # stim_probe.firstWidth = 60e-6
    # stim_probe.firstDelay = 1.0
    # stim_probe.secondLevel = 0.8
    # stim_probe.secondDelay = 10e-3
    # stim_probe.secondWidth = 60e-6
    # tab_probe = moose.Table('tab_probe')
    # tab_probe.stepMode = 3
    # tab_probe.connect('inputRequest', stim_probe, 'output')
    # stim_gate.connect('outputSrc', stim_probe, 'input')
    moose.context.setClock(0, 0.25e-3)
    moose.context.reset()
    moose.context.step(5.0)
    pylab.plot(numpy.linspace(0, 5.0, len(table_gate)), table_gate, 'rx', label='gate')
    pylab.plot(numpy.linspace(0, 5.0, len(table_gate)), tab_bg, 'b-', label='bg')
    # pylab.plot(tab_probe, 'g-', label='probe')
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    test_stim()
    


# 
# test_stim.py ends here
