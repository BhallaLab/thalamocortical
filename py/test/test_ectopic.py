#!/usr/bin/env python
# test_ectopic.py --- 
# 
# Filename: test_ectopic.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Thu May 12 11:33:31 2011 (+0530)
# Version: 
# Last-Updated: Fri May 13 15:52:05 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 66
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# test random spike generation
# 
# 

# Change log:
# 
# 
# 

# Code:

import moose
import numpy

def test_randomspike():
    count = 10
    rate = 1.0
    simdt = 1.0
    spikegens = []
    tables = []
    spike_tables = []
    model = moose.Neutral('/model')
    data = moose.Neutral('/data')
    out = open('test.out', 'w')
    for ii in range(count):
        spikegens.append(moose.RandomSpike('randspike%d' % (ii), model))
        spikegens[-1].rate = rate
        spikegens[-1].minAmp = 0.4e-9
        spikegens[-1].maxAmp = 0.4e-9
        spikegens[-1].reset = True
        spikegens[-1].resetValue = 0.0
        tables.append(moose.Table('table%d' % (ii), data))
        tables[-1].stepMode = moose.TAB_BUF
        spikegens[-1].connect('lastEvent', tables[-1], 'inputRequest')
        spike_tables.append(moose.Table('xspike%d' % (ii), data))
        spike_tables[-1].stepMode = 4
        spike_tables[-1].threshold = 0.0
        spikegens[-1].connect('state', spike_tables[-1], 'inputRequest')

        
    moose.context.setClock(0, simdt)
    moose.context.setClock(1, simdt)
    moose.context.setClock(2, simdt)
    moose.context.reset()
    for ii in range(1000):
        moose.context.step(1)
        for jj in range(count):
            out.write('%s %g\n' % (spikegens[jj].name, spikegens[jj].state))
            out.write('%s %g\n' % (spike_tables[jj].name, spike_tables[jj].output))
            jj += 1
            
    out.close()
    for child in data.children():
        table = moose.Table(child)
        table.dumpFile(table.name)
    
    
if __name__ == '__main__':
    test_randomspike()
    print 'Finished'


#  test_ectopic.py ends here
