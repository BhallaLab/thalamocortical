# synintegration.py --- 
# 
# Filename: synintegration.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Apr 24 15:30:05 2012 (+0530)
# Version: 
# Last-Updated: Tue Apr 24 15:35:48 2012 (+0530)
#           By: subha
#     Update #: 9
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This script is to set up a single spiny stellate cell and deliver
# artificial synaptic inputs to study the synaptic integration time. I
# noticed in network model that the spinystellate cells were
# responding about 40 ms after a stimulus. The purpose of this
# simulation is to find out the reason for that.
# 
# 

# Code:

from spinystellate SpinyStellate

def test_synintegration():
    sim = Simulation('synintegration')
    cell = SpinyStellate(SpinyStellate.prototype, 
                         sim.model.path + '/SpinyStellate')
    


# 
# synintegration.py ends here
