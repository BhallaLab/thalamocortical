# test_twocomp.py --- 
# 
# Filename: test_twocomp.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Tue Jun  2 23:16:38 2009 (+0530)
# Version: 
# Last-Updated: Fri Jun  5 13:58:51 2009 (+0530)
#           By: subhasis ray
#     Update #: 47
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

# Code:
import moose
import config
from compartment import MyCompartment
from simulation import Simulation

def test_twocomp(sim):
    """Pass a simulation object as parameter"""
    spine_area_mult = 2.0
    cell = moose.Cell('mycell', sim.model)
    soma = MyCompartment('soma', cell)
    soma.length = 20e-6
    soma.diameter = 2 * 7.5e-6
    soma.Em = -65e-3
    soma.initVm = -65e-3
    soma.setSpecificCm(9e-3)
    soma.setSpecificRm(5.0)
    soma.setSpecificRa(2.5)
    
    dend = MyCompartment('dend', cell)
    dend.length = 40e-6
    dend.diameter = 2 * 1.06e-6
    dend.Em = -65e-3
    dend.initVm = -65e-3
    dend.setSpecificCm(9e-3 * spine_area_mult)
    dend.setSpecificRm(5.0 / spine_area_mult)
    dend.setSpecificRa(2.5)
    soma.traubConnect(dend)
    vm_table = dend.insertRecorder("Vm1", "Vm", sim.data)
    dend.insertPulseGen("pulsegen1", sim.model, firstLevel=3e-10, firstDelay=20e-3, firstWidth=100e-3)
#     sim.schedule()
#     sim.run(100e-3)
#     sim.dump_data("data")
    return cell
# print 'soma:', 'Ra =', soma.Ra, 'Rm =', soma.Rm, 'Cm =', soma.Cm, 'Em =', soma.Em, 'initVm =', soma.initVm
# print 'dend:', 'Ra =', dend.Ra, 'Rm =', dend.Rm, 'Cm =', dend.Cm, 'Em =', dend.Em, 'initVm =', dend.initVm
# 
# test_twocomp.py ends here
