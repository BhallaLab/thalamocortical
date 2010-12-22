# test_stim.py --- 
# 
# Filename: test_stim.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Wed Dec 22 17:29:05 2010 (+0530)
# Version: 
# Last-Updated: Wed Dec 22 22:04:21 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 95
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

import pylab
import numpy
import moose

class StimTest:
    def __init__(self):
        stim_container = moose.Neutral('stim')
        self.setup_stimulus(stim_container, 0.5, 0.1, 0.01, 0.001, 0.005)
        self.tables = {}
        tab_cnt = 0
        self.data = moose.Neutral('data')
        attrs = ['firstLevel', 'firstDelay', 'firstWidth', 'secondLevel', 'secondDelay', 'secondWidth', 'trigMode']
        for stim_id in self.stim_container.children():
            stim = moose.PulseGen(stim_id)
            print stim.path
            for attr in attrs:
                print attr, eval('stim.' + attr)
            tab = moose.Table(stim.name, self.data)
            tab.stepMode = 3
            tab.connect('inputRequest', stim, 'output')
            self.tables[stim.name] = tab
            tab_cnt += 1

        moose.context.setClock(0, 1e-4)
        moose.context.setClock(1, 1e-4)
        moose.context.setClock(2, 1e-4)
        moose.context.setClock(3, 1e-4)
        moose.context.reset()
        moose.context.step(1.0)
        cnt = 1
        x_data = numpy.linspace(0, 1.0, len(self.tables['root_gate']))
        pylab.subplot(211)
        pylab.plot(x_data, numpy.array(self.tables['stim_background']), label='background')
        pylab.plot(x_data, numpy.array(self.tables['stim_probe']), 'r-.', label='probe')
        pylab.legend(loc='upper left')
        pylab.subplot(212)
        pylab.plot(x_data, numpy.array(self.tables['root_gate']), label='onset-gate')
        pylab.plot(x_data, numpy.array(self.tables['stim_gate']), 'g-.', label='stimulus gate')
        pylab.legend(loc='upper left')                               
        pylab.show()

    def setup_stimulus(self, stim_container, stim_onset, stim_interval, bg_delay, pulse_width, isi, level=5e-12, bg_count=10):
        """Setup the stimulus protocol.

        The protocol is as follows:

        Let the system stabilize for stim_onset seconds.

        Then turn the stim_gate on: which gates the triggers.

        The bg_trigger will trigger the background pulse
        generator. probe_trigger will trigger the probe pulse
        generator.

        stim_container -- container object for stimulating electrodes

        celltype -- type of cells we are looking at

        stim_onset -- when we consider the system stabilized and start
        stimulus

        stim_interval -- interval between two stimulus sessions

        bg_delay -- start time for background stimulus after
        stim_onset, this is also the time between successive
        applications of background.

        pulse_width -- width of background pulses

        isi -- if paired pulse, then the interval between the
        two pulses (beginning of second - beginning of first), 0 for
        single pulse.

        level -- current injection value

        bg_count -- number of cells stimulated by background pulse.

        """
        if  isinstance(stim_container, str):
            self.stim_container = moose.Neutral(stim_container)
        elif isinstance(stim_container, moose.Neutral):
            self.stim_container = stim_container
        else:
            print isinstance(stim_container, moose.Neutral)
            raise Exception('Stimulus container must be a string or a Neutral object: got %s', stim_container.__class__.__name__)
        self.root_gate = moose.PulseGen('root_gate', self.stim_container)
        self.root_gate.trigMode = 0 # Free running
        self.root_gate.firstDelay = stim_onset
        self.root_gate.firstWidth = 1e9 # Keep it on forever
        self.root_gate.firstLevel = 1.0                
        self.stim_gate = moose.PulseGen('stim_gate', self.stim_container)
        self.stim_gate.firstDelay = stim_interval
        self.stim_gate.firstLevel = 1.0
        self.stim_gate.firstWidth = (bg_delay + pulse_width + isi) * 2
        self.stim_gate.trigMode = 2
        self.stim_bg = moose.PulseGen('stim_background', self.stim_container)
        self.stim_bg.firstLevel = level
        self.stim_bg.secondLevel = level
        self.stim_bg.firstDelay = bg_delay
        self.stim_bg.firstWidth = pulse_width
        self.stim_bg.secondDelay = isi
        self.stim_bg.secondWidth = pulse_width
        
        self.stim_bg.trigMode = 2
        self.stim_probe = moose.PulseGen('stim_probe', self.stim_container)
        self.stim_probe.firstLevel = level
        self.stim_probe.secondLevel = level
        self.stim_probe.firstDelay = 2 * bg_delay + pulse_width + isi
        self.stim_probe.secondDelay = isi
        self.stim_probe.firstWidth = pulse_width
        self.stim_probe.secondWidth = pulse_width            
        self.stim_probe.trigMode = 2
        self.root_gate.connect('outputSrc', self.stim_gate, 'input')
        self.stim_gate.connect('outputSrc', self.stim_bg, 'input')
        self.stim_gate.connect('outputSrc', self.stim_probe, 'input')


if __name__ == '__main__':
    a = StimTest()
            

# 
# test_stim.py ends here
