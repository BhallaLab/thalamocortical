# test_pulse.py --- 
# 
# Filename: test_pulse.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Mon Dec 20 11:42:39 2010 (+0530)
# Version: 
# Last-Updated: Mon Dec 20 16:19:16 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 36
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
# How to do the stimulus?
# 1. Gated mode of background:
#  stabilization gate  (on from 0.5 s) -> background_trigger
# background_trigger (frequency = f_t) -> all_backgrounds (2-pulses per on state)
# probe_trigger (frequency= 2*f_t) -> probe (2-pulses per on state)

#
# I have a strong feeling that to void any pulses while the system is
# stabilizing, we should allow an additional initial delay (may be
# call it baseDelay?)

trigger = moose.PulseGen('trigger')
trigger.trigMode = 0 # Free run
trigger.firstLevel = 1.0
trigger.firstDelay = 50.0
trigger.firstWidth = 10.0

pulse = moose.PulseGen('pulse')
pulse.trigMode = 2 # external gate
pulse.firstLevel = -10.0
pulse.firstWidth = 2.0
pulse.firstDelay = 2.0

trigger.connect('outputSrc', pulse, 'input')

trigtab = moose.Table('trigtab')
trigtab.stepMode = 3
trigtab.connect('inputRequest', trigger, 'output')

pulsetab = moose.Table('pulsetab')
pulsetab.stepMode = 3
pulsetab.connect('inputRequest', pulse, 'output')

moose.context.reset()
moose.context.step(1000.0)

pylab.plot(numpy.array(trigtab), label='Trigger')
pylab.plot(numpy.array(pulsetab), label='Pulse')
pylab.legend()
pylab.show()




# 
# test_pulse.py ends here
