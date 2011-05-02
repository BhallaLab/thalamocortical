# synapse.py --- 
# 
# Filename: synapse.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Feb 25 15:22:11 2010 (+0530)
# Version: 
# Last-Updated: Mon May  2 11:24:13 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 183
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This file is for data on synapses of various
# kinds. Ideally it should be replaced by a datafile (something in
# netCDF/HDF5?). But given the small size of this data, I don't see
# much savings in that.
# 
#

# Change log:
# 
#  2011-04-29 12:32:24 (+0530) - adding Pr and delay 
# 
# 

# Code:

from collections import defaultdict
import numpy

E_GABA = {
    'nRT': -75e-3,
    'SupBasket': -75e-3,
    'SupAxoaxonic': -75e-3,
    'SupLTS': -75e-3,
    # 'SupPyrFRB': -81e-3,
    # 'SupPyrRS': -81e-3,
    'DeepBasket': -75e-3,
    'DeepAxoaxonic': -75e-3,
    'DeepLTS': -75e-3,
    'nRT': -75e-3,
    # 'SpinyStellate': -75e-3,
    # 'NontuftedRS': -75e-3,
    # 'TuftedIB': -75e-3,
    # 'TuftedRS': -75e-3
}
TAU_AMPA = defaultdict(dict)
TAU_NMDA = defaultdict(dict)
TAU_GABA = defaultdict(dict)
TAU_GABA_FAST = defaultdict(dict)
TAU_GABA_SLOW = defaultdict(dict)
G_AMPA = defaultdict(dict)
G_NMDA = defaultdict(dict)
G_GABA = defaultdict(dict)



# Synaptic conductance time constants. 
TAU_AMPA['SupPyrRS']['SupPyrRS'] = 2.e0 
TAU_NMDA['SupPyrRS']['SupPyrRS'] = 130.5e0 
TAU_AMPA['SupPyrRS']['SupPyrFRB'] = 2.e0 
TAU_NMDA['SupPyrRS']['SupPyrFRB'] = 130.e0  
TAU_AMPA['SupPyrRS']['SupBasket'] = .8e0   
TAU_NMDA['SupPyrRS']['SupBasket'] = 100.e0 
TAU_AMPA['SupPyrRS']['SupAxoaxonic'] = .8e0  
TAU_NMDA['SupPyrRS']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['SupPyrRS']['SupLTS'] = 1.e0  
TAU_NMDA['SupPyrRS']['SupLTS'] = 100.e0 
TAU_AMPA['SupPyrRS']['SpinyStellate'] = 2.e0   
TAU_NMDA['SupPyrRS']['SpinyStellate'] = 130.e0 
TAU_AMPA['SupPyrRS']['TuftedIB'] = 2.e0   
TAU_NMDA['SupPyrRS']['TuftedIB'] = 130.e0 
TAU_AMPA['SupPyrRS']['TuftedRS'] = 2.e0   
TAU_NMDA['SupPyrRS']['TuftedRS'] = 130.e0 
TAU_AMPA['SupPyrRS']['DeepBasket'] = .8e0   
TAU_NMDA['SupPyrRS']['DeepBasket'] = 100.e0 
TAU_AMPA['SupPyrRS']['DeepAxoaxonic'] = .8e0   
TAU_NMDA['SupPyrRS']['DeepAxoaxonic'] = 100.e0 
TAU_AMPA['SupPyrRS']['DeepLTS'] = 1.e0   
TAU_NMDA['SupPyrRS']['DeepLTS'] = 100.e0 
TAU_AMPA['SupPyrRS']['NontuftedRS'] = 2.e0   
TAU_NMDA['SupPyrRS']['NontuftedRS'] = 130.e0 

TAU_AMPA['SupPyrFRB']['SupPyrRS'] = 2.e0   
TAU_NMDA['SupPyrFRB']['SupPyrRS'] = 130.e0  
TAU_AMPA['SupPyrFRB']['SupPyrFRB'] = 2.e0   
TAU_NMDA['SupPyrFRB']['SupPyrFRB'] = 130.e0 
TAU_AMPA['SupPyrFRB']['SupBasket'] = .8e0   
TAU_NMDA['SupPyrFRB']['SupBasket'] = 100.e0 
TAU_AMPA['SupPyrFRB']['SupAxoaxonic'] = .8e0  
TAU_NMDA['SupPyrFRB']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['SupPyrFRB']['SupLTS'] = 1.e0  
TAU_NMDA['SupPyrFRB']['SupLTS'] = 100.e0 
TAU_AMPA['SupPyrFRB']['SpinyStellate'] = 2.e0   
TAU_NMDA['SupPyrFRB']['SpinyStellate'] = 130.e0 
TAU_AMPA['SupPyrFRB']['TuftedIB'] = 2.e0   
TAU_NMDA['SupPyrFRB']['TuftedIB'] = 130.e0 
TAU_AMPA['SupPyrFRB']['TuftedRS'] = 2.e0   
TAU_NMDA['SupPyrFRB']['TuftedRS'] = 130.e0
TAU_AMPA['SupPyrFRB']['DeepBasket'] = .8e0   
TAU_NMDA['SupPyrFRB']['DeepBasket'] = 100.e0 
TAU_AMPA['SupPyrFRB']['DeepAxoaxonic'] = .8e0   
TAU_NMDA['SupPyrFRB']['DeepAxoaxonic'] = 100.e0 
TAU_AMPA['SupPyrFRB']['DeepLTS'] = 1.e0   
TAU_NMDA['SupPyrFRB']['DeepLTS'] = 100.e0 
TAU_AMPA['SupPyrFRB']['NontuftedRS'] = 2.e0   
TAU_NMDA['SupPyrFRB']['NontuftedRS'] = 130.e0 

TAU_GABA['SupBasket']['SupPyrRS'] = 6.e0  
TAU_GABA['SupBasket']['SupPyrFRB'] = 6.e0  
TAU_GABA['SupBasket']['SupBasket'] = 3.e0  
TAU_GABA['SupBasket']['SupAxoaxonic'] = 3.e0  
TAU_GABA['SupBasket']['SupLTS'] = 3.e0  
TAU_GABA['SupBasket']['SpinyStellate'] = 6.e0  

TAU_GABA['SupAxoaxonic']['SupPyrRS'] = 6.e0  
TAU_GABA['SupAxoaxonic']['SupPyrFRB'] = 6.e0  
TAU_GABA['SupAxoaxonic']['SpinyStellate'] = 6.e0  
TAU_GABA['SupAxoaxonic']['TuftedIB'] = 6.e0  
TAU_GABA['SupAxoaxonic']['TuftedRS'] = 6.e0  
TAU_GABA['SupAxoaxonic']['NontuftedRS'] = 6.e0  

TAU_GABA['SupLTS']['SupPyrRS'] = 20.e0 
TAU_GABA['SupLTS']['SupPyrFRB'] = 20.e0 
TAU_GABA['SupLTS']['SupBasket'] = 20.e0 
TAU_GABA['SupLTS']['SupAxoaxonic'] = 20.e0 
TAU_GABA['SupLTS']['SupLTS'] = 20.e0 
TAU_GABA['SupLTS']['SpinyStellate'] = 20.e0 
TAU_GABA['SupLTS']['TuftedIB'] = 20.e0 
TAU_GABA['SupLTS']['TuftedRS'] = 20.e0 
TAU_GABA['SupLTS']['DeepBasket'] = 20.e0 
TAU_GABA['SupLTS']['DeepAxoaxonic'] = 20.e0 
TAU_GABA['SupLTS']['DeepLTS'] = 20.e0 
TAU_GABA['SupLTS']['NontuftedRS'] = 20.e0  

TAU_AMPA['SpinyStellate']['SupPyrRS'] = 2.e0  
TAU_NMDA['SpinyStellate']['SupPyrRS'] = 130.e0 
TAU_AMPA['SpinyStellate']['SupPyrFRB'] = 2.e0  
TAU_NMDA['SpinyStellate']['SupPyrFRB'] = 130.e0 
TAU_AMPA['SpinyStellate']['SupBasket'] = .8e0  
TAU_NMDA['SpinyStellate']['SupBasket'] = 100.e0
TAU_AMPA['SpinyStellate']['SupAxoaxonic'] = .8e0  
TAU_NMDA['SpinyStellate']['SupAxoaxonic'] = 100.e0
TAU_AMPA['SpinyStellate']['SupLTS'] = 1.e0  
TAU_NMDA['SpinyStellate']['SupLTS'] = 100.e0
TAU_AMPA['SpinyStellate']['SpinyStellate'] = 2.e0  
TAU_NMDA['SpinyStellate']['SpinyStellate'] = 130.e0 
TAU_AMPA['SpinyStellate']['TuftedIB'] = 2.e0  
TAU_NMDA['SpinyStellate']['TuftedIB'] = 130.e0 
TAU_AMPA['SpinyStellate']['TuftedRS'] = 2.e0  
TAU_NMDA['SpinyStellate']['TuftedRS'] = 130.e0
TAU_AMPA['SpinyStellate']['DeepBasket'] = .8e0  
TAU_NMDA['SpinyStellate']['DeepBasket'] = 100.e0
TAU_AMPA['SpinyStellate']['DeepAxoaxonic'] = .8e0  
TAU_NMDA['SpinyStellate']['DeepAxoaxonic'] = 100.e0
TAU_AMPA['SpinyStellate']['DeepLTS'] = 1.e0  
TAU_NMDA['SpinyStellate']['DeepLTS'] = 100.e0
TAU_AMPA['SpinyStellate']['NontuftedRS'] = 2.e0  
TAU_NMDA['SpinyStellate']['NontuftedRS'] = 130.e0

TAU_AMPA['TuftedIB']['SupPyrRS'] = 2.e0 
TAU_NMDA['TuftedIB']['SupPyrRS'] = 130.e0
TAU_AMPA['TuftedIB']['SupPyrFRB'] = 2.e0 
TAU_NMDA['TuftedIB']['SupPyrFRB'] = 130.e0
TAU_AMPA['TuftedIB']['SupBasket'] = .8e0  
TAU_NMDA['TuftedIB']['SupBasket'] = 100.e0 
TAU_AMPA['TuftedIB']['SupAxoaxonic'] = .8e0  
TAU_NMDA['TuftedIB']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['TuftedIB']['SupLTS'] = 1.e0  
TAU_NMDA['TuftedIB']['SupLTS'] = 100.e0 
TAU_AMPA['TuftedIB']['SpinyStellate'] = 2.e0   
TAU_NMDA['TuftedIB']['SpinyStellate'] = 130.e0 
TAU_AMPA['TuftedIB']['TuftedIB'] = 2.e0  
TAU_NMDA['TuftedIB']['TuftedIB'] = 130.e0 
TAU_AMPA['TuftedIB']['TuftedRS'] = 2.0e0 
TAU_NMDA['TuftedIB']['TuftedRS'] = 130.e0 
TAU_AMPA['TuftedIB']['DeepBasket'] = .8e0  
TAU_NMDA['TuftedIB']['DeepBasket'] = 100.e0 
TAU_AMPA['TuftedIB']['DeepAxoaxonic'] = .8e0  
TAU_NMDA['TuftedIB']['DeepAxoaxonic'] = 100.e0 
TAU_AMPA['TuftedIB']['DeepLTS'] = 1.e0  
TAU_NMDA['TuftedIB']['DeepLTS'] = 100.e0 
TAU_AMPA['TuftedIB']['NontuftedRS'] = 2.0e0 
TAU_NMDA['TuftedIB']['NontuftedRS'] = 130.e0 

TAU_AMPA['TuftedRS']['SupPyrRS'] = 2.e0 
TAU_NMDA['TuftedRS']['SupPyrRS'] = 130.e0
TAU_AMPA['TuftedRS']['SupPyrFRB'] = 2.e0 
TAU_NMDA['TuftedRS']['SupPyrFRB'] = 130.e0
TAU_AMPA['TuftedRS']['SupBasket'] = .8e0  
TAU_NMDA['TuftedRS']['SupBasket'] = 100.e0 
TAU_AMPA['TuftedRS']['SupAxoaxonic'] = .8e0  
TAU_NMDA['TuftedRS']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['TuftedRS']['SupLTS'] = 1.e0  
TAU_NMDA['TuftedRS']['SupLTS'] = 100.e0 
TAU_AMPA['TuftedRS']['SpinyStellate'] = 2.e0  
TAU_NMDA['TuftedRS']['SpinyStellate'] = 130.e0 
TAU_AMPA['TuftedRS']['TuftedIB'] = 2.e0  
TAU_NMDA['TuftedRS']['TuftedIB'] = 130.e0 
TAU_AMPA['TuftedRS']['TuftedRS'] = 2.e0  
TAU_NMDA['TuftedRS']['TuftedRS'] = 130.e0 
TAU_AMPA['TuftedRS']['DeepBasket'] = .8e0  
TAU_NMDA['TuftedRS']['DeepBasket'] = 100.e0 
TAU_AMPA['TuftedRS']['DeepAxoaxonic'] = .8e0  
TAU_NMDA['TuftedRS']['DeepAxoaxonic'] = 100.e0 
TAU_AMPA['TuftedRS']['DeepLTS'] = 1.e0   
TAU_NMDA['TuftedRS']['DeepLTS'] = 100.e0 
TAU_AMPA['TuftedRS']['NontuftedRS'] = 2.e0  
TAU_NMDA['TuftedRS']['NontuftedRS'] = 130.e0 

TAU_GABA['DeepBasket']['SpinyStellate'] = 6.e0  
TAU_GABA['DeepBasket']['TuftedIB'] = 6.e0  
TAU_GABA['DeepBasket']['TuftedRS'] = 6.e0  
TAU_GABA['DeepBasket']['DeepBasket'] = 3.e0  
TAU_GABA['DeepBasket']['DeepAxoaxonic'] = 3.e0  
TAU_GABA['DeepBasket']['DeepLTS'] = 3.e0  
TAU_GABA['DeepBasket']['NontuftedRS'] = 6.e0  

TAU_GABA['DeepAxoaxonic']['SupPyrRS'] = 6.e0  
TAU_GABA['DeepAxoaxonic']['SupPyrFRB'] = 6.e0  
TAU_GABA['DeepAxoaxonic']['SpinyStellate'] = 6.e0  
TAU_GABA['DeepAxoaxonic']['TuftedIB'] = 6.e0  
TAU_GABA['DeepAxoaxonic']['TuftedRS'] = 6.e0  
TAU_GABA['DeepAxoaxonic']['NontuftedRS'] = 6.e0  

TAU_GABA['DeepLTS']['SupPyrRS'] = 20.e0 
TAU_GABA['DeepLTS']['SupPyrFRB'] = 20.e0 
TAU_GABA['DeepLTS']['SupBasket'] = 20.e0 
TAU_GABA['DeepLTS']['SupAxoaxonic'] = 20.e0 
TAU_GABA['DeepLTS']['SupLTS'] = 20.e0 
TAU_GABA['DeepLTS']['SpinyStellate'] = 20.e0 
TAU_GABA['DeepLTS']['TuftedIB'] = 20.e0 
TAU_GABA['DeepLTS']['TuftedRS'] = 20.e0 
TAU_GABA['DeepLTS']['DeepBasket'] = 20.e0 
TAU_GABA['DeepLTS']['DeepAxoaxonic'] = 20.e0 
TAU_GABA['DeepLTS']['DeepLTS'] = 20.e0 
TAU_GABA['DeepLTS']['NontuftedRS'] = 20.e0 

TAU_AMPA['TCR']['SupPyrRS'] = 2.e0  
TAU_NMDA['TCR']['SupPyrRS'] = 130.e0
TAU_AMPA['TCR']['SupPyrFRB'] = 2.e0  
TAU_NMDA['TCR']['SupPyrFRB'] = 130.e0
TAU_AMPA['TCR']['SupBasket'] = 1.e0  
TAU_NMDA['TCR']['SupBasket'] = 100.e0
TAU_AMPA['TCR']['SupAxoaxonic'] = 1.e0  
TAU_NMDA['TCR']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['TCR']['SpinyStellate'] = 2.0e0 
TAU_NMDA['TCR']['SpinyStellate'] = 130.e0
TAU_AMPA['TCR']['TuftedIB'] = 2.e0  
TAU_NMDA['TCR']['TuftedIB'] = 130.e0
TAU_AMPA['TCR']['TuftedRS'] = 2.e0  
TAU_NMDA['TCR']['TuftedRS'] = 130.e0
TAU_AMPA['TCR']['DeepBasket'] = 1.e0  
TAU_NMDA['TCR']['DeepBasket'] = 100.e0
TAU_AMPA['TCR']['DeepAxoaxonic'] = 1.e0  
TAU_NMDA['TCR']['DeepAxoaxonic'] = 100.e0
TAU_AMPA['TCR']['nRT'] = 2.0e0      
TAU_NMDA['TCR']['nRT'] = 150.e0
TAU_AMPA['TCR']['NontuftedRS'] = 2.0e0     
TAU_NMDA['TCR']['NontuftedRS'] = 130.e0

#         TAU_GABA_FAST['nRT']['TCR'] = 10.e0 
#         TAU_GABA_SLOW['nRT']['TCR'] = 30.e0 
#         TAU_GABA_FAST['nRT']['nRT'] = 18.e0 
#         TAU_GABA_SLOW['nRT']['nRT'] = 89.e0 
# See notebook entry of 17 Feb. 2004.
# Speed these up per Huntsman & Huguenard (2000)
TAU_GABA_FAST['nRT']['TCR'] = 3.30e0 
TAU_GABA_SLOW['nRT']['TCR'] = 10.e0 
TAU_GABA_FAST['nRT']['nRT'] = 9.e0 
TAU_GABA_SLOW['nRT']['nRT'] = 44.5e0 

TAU_AMPA['NontuftedRS']['SupPyrRS'] = 2.e0  
TAU_NMDA['NontuftedRS']['SupPyrRS'] = 130.e0
TAU_AMPA['NontuftedRS']['SupPyrFRB'] = 2.e0  
TAU_NMDA['NontuftedRS']['SupPyrFRB'] = 130.e0
TAU_AMPA['NontuftedRS']['SupBasket'] = .8e0  
TAU_NMDA['NontuftedRS']['SupBasket'] = 100.e0
TAU_AMPA['NontuftedRS']['SupAxoaxonic'] = .8e0  
TAU_NMDA['NontuftedRS']['SupAxoaxonic'] = 100.e0 
TAU_AMPA['NontuftedRS']['SupLTS'] = 1.0e0 
TAU_NMDA['NontuftedRS']['SupLTS'] = 100.e0
TAU_AMPA['NontuftedRS']['SpinyStellate'] = 2.e0  
TAU_NMDA['NontuftedRS']['SpinyStellate'] = 130.e0
TAU_AMPA['NontuftedRS']['TuftedIB'] = 2.e0  
TAU_NMDA['NontuftedRS']['TuftedIB'] = 130.e0
TAU_AMPA['NontuftedRS']['TuftedRS'] = 2.e0  
TAU_NMDA['NontuftedRS']['TuftedRS'] = 130.e0
TAU_AMPA['NontuftedRS']['DeepBasket'] = .8e0  
TAU_NMDA['NontuftedRS']['DeepBasket'] = 100.e0
TAU_AMPA['NontuftedRS']['DeepAxoaxonic'] = .8e0   
TAU_NMDA['NontuftedRS']['DeepAxoaxonic'] = 100.e0

TAU_AMPA['NontuftedRS']['DeepLTS'] = 1.e0  
TAU_NMDA['NontuftedRS']['DeepLTS'] = 100.e0
TAU_AMPA['NontuftedRS']['TCR'] = 2.e0  
TAU_NMDA['NontuftedRS']['TCR'] = 130.e0 
TAU_AMPA['NontuftedRS']['nRT'] = 2.0e0 
TAU_NMDA['NontuftedRS']['nRT'] = 100.e0 
TAU_AMPA['NontuftedRS']['NontuftedRS'] = 2.e0  
TAU_NMDA['NontuftedRS']['NontuftedRS'] = 130.e0 
# End definition of synaptic time constants

# Synaptic conductance scaling factors.
# The values are taken directly from the NEURON file: groucho.hoc
# The units are updated to SI programmatically later ( uS -> S)
G_AMPA['SupPyrRS']['SupPyrRS'] = 0.25e-3
G_NMDA['SupPyrRS']['SupPyrRS'] = 0.025e-3
G_AMPA['SupPyrRS']['SupPyrFRB'] = 0.25e-3
G_NMDA['SupPyrRS']['SupPyrFRB'] = 0.025e-3
G_AMPA['SupPyrRS']['SupBasket'] = 3.00e-3
G_NMDA['SupPyrRS']['SupBasket'] = 0.15e-3
G_AMPA['SupPyrRS']['SupAxoaxonic'] = 3.0e-3
G_NMDA['SupPyrRS']['SupAxoaxonic'] = 0.15e-3
G_AMPA['SupPyrRS']['SupLTS'] = 2.0e-3
G_NMDA['SupPyrRS']['SupLTS'] = 0.15e-3
G_AMPA['SupPyrRS']['SpinyStellate'] = 0.10e-3
G_NMDA['SupPyrRS']['SpinyStellate'] = 0.01e-3
G_AMPA['SupPyrRS']['TuftedIB'] = 0.10e-3
G_NMDA['SupPyrRS']['TuftedIB'] = 0.01e-3
G_AMPA['SupPyrRS']['TuftedRS'] = 0.10e-3
G_NMDA['SupPyrRS']['TuftedRS'] = 0.01e-3
G_AMPA['SupPyrRS']['DeepBasket'] = 1.00e-3
G_NMDA['SupPyrRS']['DeepBasket'] = 0.10e-3
G_AMPA['SupPyrRS']['DeepAxoaxonic'] = 1.00e-3
G_NMDA['SupPyrRS']['DeepAxoaxonic'] = 0.10e-3
G_AMPA['SupPyrRS']['DeepLTS'] = 1.00e-3
G_NMDA['SupPyrRS']['DeepLTS'] = 0.15e-3
G_AMPA['SupPyrRS']['NontuftedRS'] = 0.50e-3
G_NMDA['SupPyrRS']['NontuftedRS'] = 0.05e-3

G_AMPA['SupPyrFRB']['SupPyrRS'] = 0.25e-3
G_NMDA['SupPyrFRB']['SupPyrRS'] = 0.025e-3
G_AMPA['SupPyrFRB']['SupPyrFRB'] = 0.25e-3
G_NMDA['SupPyrFRB']['SupPyrFRB'] = .025e-3
G_AMPA['SupPyrFRB']['SupBasket'] = 3.00e-3
G_NMDA['SupPyrFRB']['SupBasket'] = 0.10e-3
G_AMPA['SupPyrFRB']['SupAxoaxonic'] = 3.0e-3
G_NMDA['SupPyrFRB']['SupAxoaxonic'] = 0.10e-3
G_AMPA['SupPyrFRB']['SupLTS'] = 2.0e-3
G_NMDA['SupPyrFRB']['SupLTS'] = 0.10e-3
G_AMPA['SupPyrFRB']['SpinyStellate'] = 0.10e-3
G_NMDA['SupPyrFRB']['SpinyStellate'] = 0.01e-3
G_AMPA['SupPyrFRB']['TuftedIB'] = 0.10e-3
G_NMDA['SupPyrFRB']['TuftedIB'] = 0.01e-3
G_AMPA['SupPyrFRB']['TuftedRS'] = 0.10e-3
G_NMDA['SupPyrFRB']['TuftedRS'] = 0.01e-3
G_AMPA['SupPyrFRB']['DeepBasket'] = 1.00e-3
G_NMDA['SupPyrFRB']['DeepBasket'] = 0.10e-3
G_AMPA['SupPyrFRB']['DeepAxoaxonic'] = 1.00e-3
G_NMDA['SupPyrFRB']['DeepAxoaxonic'] = 0.10e-3
G_AMPA['SupPyrFRB']['DeepLTS'] = 1.00e-3
G_NMDA['SupPyrFRB']['DeepLTS'] = 0.10e-3
G_AMPA['SupPyrFRB']['NontuftedRS'] = 0.50e-3
G_NMDA['SupPyrFRB']['NontuftedRS'] = 0.05e-3

G_GABA['SupBasket']['SupPyrRS'] = 1.2e-3
G_GABA['SupBasket']['SupPyrFRB'] = 1.2e-3
G_GABA['SupBasket']['SupBasket'] = 0.2e-3
G_GABA['SupBasket']['SupAxoaxonic'] = 0.2e-3
G_GABA['SupBasket']['SupLTS'] = 0.5e-3
#      G_GABA['SupBasket']['SpinyStellate'] = 0.7e-3
G_GABA['SupBasket']['SpinyStellate'] = 0.1e-3 # if main inhib. to spinstell from deep int.

G_GABA['SupAxoaxonic']['SupPyrRS'] = 1.2e-3
G_GABA['SupAxoaxonic']['SupPyrFRB'] = 1.2e-3
#      G_GABA['SupAxoaxonic']['SpinyStellate'] = 1.0e-3
G_GABA['SupAxoaxonic']['SpinyStellate'] = 0.1e-3 # if main inhib. to spinstell from deep int.
G_GABA['SupAxoaxonic']['TuftedIB'] = 1.0e-3
G_GABA['SupAxoaxonic']['TuftedRS'] = 1.0e-3
G_GABA['SupAxoaxonic']['NontuftedRS'] = 1.0e-3


G_GABA['SupLTS']['SupPyrRS'] = .01e-3
G_GABA['SupLTS']['SupPyrFRB'] = .01e-3
G_GABA['SupLTS']['SupBasket'] = .01e-3
G_GABA['SupLTS']['SupAxoaxonic'] = .01e-3
G_GABA['SupLTS']['SupLTS'] = .05e-3
G_GABA['SupLTS']['SpinyStellate'] = .01e-3
G_GABA['SupLTS']['TuftedIB'] = .02e-3
G_GABA['SupLTS']['TuftedRS'] = .02e-3
G_GABA['SupLTS']['DeepBasket'] = .01e-3
G_GABA['SupLTS']['DeepAxoaxonic'] = .01e-3
G_GABA['SupLTS']['DeepLTS'] = .05e-3
G_GABA['SupLTS']['NontuftedRS'] = .01e-3

G_AMPA['SpinyStellate']['SupPyrRS'] = 1.0e-3
G_NMDA['SpinyStellate']['SupPyrRS'] = 0.1e-3
G_AMPA['SpinyStellate']['SupPyrFRB'] = 1.0e-3
G_NMDA['SpinyStellate']['SupPyrFRB'] = 0.1e-3
G_AMPA['SpinyStellate']['SupBasket'] = 1.0e-3
G_NMDA['SpinyStellate']['SupBasket'] = .15e-3
G_AMPA['SpinyStellate']['SupAxoaxonic'] = 1.0e-3
G_NMDA['SpinyStellate']['SupAxoaxonic'] = .15e-3
G_AMPA['SpinyStellate']['SupLTS'] = 1.0e-3
G_NMDA['SpinyStellate']['SupLTS'] = .15e-3
G_AMPA['SpinyStellate']['SpinyStellate'] = 1.0e-3
G_NMDA['SpinyStellate']['SpinyStellate'] = 0.1e-3
G_AMPA['SpinyStellate']['TuftedIB'] = 1.0e-3
G_NMDA['SpinyStellate']['TuftedIB'] = 0.1e-3
G_AMPA['SpinyStellate']['TuftedRS'] = 1.0e-3
G_NMDA['SpinyStellate']['TuftedRS'] = 0.1e-3
G_AMPA['SpinyStellate']['DeepBasket'] = 1.0e-3
G_NMDA['SpinyStellate']['DeepBasket'] = .15e-3
G_AMPA['SpinyStellate']['DeepAxoaxonic'] = 1.0e-3
G_NMDA['SpinyStellate']['DeepAxoaxonic'] = .15e-3
G_AMPA['SpinyStellate']['DeepLTS'] = 1.0e-3
G_NMDA['SpinyStellate']['DeepLTS'] = .15e-3
G_AMPA['SpinyStellate']['NontuftedRS'] = 1.0e-3
G_NMDA['SpinyStellate']['NontuftedRS'] = 0.1e-3

G_AMPA['TuftedIB']['SupPyrRS'] = 0.5e-3
G_NMDA['TuftedIB']['SupPyrRS'] = 0.05e-3
G_AMPA['TuftedIB']['SupPyrFRB'] = 0.5e-3
G_NMDA['TuftedIB']['SupPyrFRB'] = 0.05e-3
G_AMPA['TuftedIB']['SupBasket'] = 1.0e-3
G_NMDA['TuftedIB']['SupBasket'] = 0.15e-3
G_AMPA['TuftedIB']['SupAxoaxonic'] = 1.0e-3
G_NMDA['TuftedIB']['SupAxoaxonic'] = 0.15e-3
G_AMPA['TuftedIB']['SupLTS'] = 1.0e-3
G_NMDA['TuftedIB']['SupLTS'] = 0.15e-3
G_AMPA['TuftedIB']['SpinyStellate'] = 0.5e-3
G_NMDA['TuftedIB']['SpinyStellate'] = 0.05e-3
G_AMPA['TuftedIB']['TuftedIB'] = 2.0e-3
G_NMDA['TuftedIB']['TuftedIB'] = 0.20e-3
G_AMPA['TuftedIB']['TuftedRS'] = 2.0e-3
G_NMDA['TuftedIB']['TuftedRS'] = 0.20e-3
G_AMPA['TuftedIB']['DeepBasket'] = 3.0e-3
G_NMDA['TuftedIB']['DeepBasket'] = 0.15e-3
G_AMPA['TuftedIB']['DeepAxoaxonic'] = 3.0e-3
G_NMDA['TuftedIB']['DeepAxoaxonic'] = 0.15e-3
G_AMPA['TuftedIB']['DeepLTS'] = 2.0e-3
G_NMDA['TuftedIB']['DeepLTS'] = 0.15e-3
G_AMPA['TuftedIB']['NontuftedRS'] = 2.0e-3
G_NMDA['TuftedIB']['NontuftedRS'] = 0.20e-3

G_AMPA['TuftedRS']['SupPyrRS'] = 0.5e-3
G_NMDA['TuftedRS']['SupPyrRS'] = 0.05e-3
G_AMPA['TuftedRS']['SupPyrFRB'] = 0.5e-3
G_NMDA['TuftedRS']['SupPyrFRB'] = 0.05e-3
G_AMPA['TuftedRS']['SupBasket'] = 1.0e-3
G_NMDA['TuftedRS']['SupBasket'] = 0.15e-3
G_AMPA['TuftedRS']['SupAxoaxonic'] = 1.0e-3
G_NMDA['TuftedRS']['SupAxoaxonic'] = 0.15e-3
G_AMPA['TuftedRS']['SupLTS'] = 1.0e-3
G_NMDA['TuftedRS']['SupLTS'] = 0.15e-3
G_AMPA['TuftedRS']['SpinyStellate'] = 0.5e-3
G_NMDA['TuftedRS']['SpinyStellate'] = 0.05e-3
G_AMPA['TuftedRS']['TuftedIB'] = 1.0e-3
G_NMDA['TuftedRS']['TuftedIB'] = 0.10e-3
G_AMPA['TuftedRS']['TuftedRS'] = 1.0e-3
G_NMDA['TuftedRS']['TuftedRS'] = 0.10e-3
G_AMPA['TuftedRS']['DeepBasket'] = 3.0e-3
G_NMDA['TuftedRS']['DeepBasket'] = 0.10e-3
G_AMPA['TuftedRS']['DeepAxoaxonic'] = 3.0e-3
G_NMDA['TuftedRS']['DeepAxoaxonic'] = 0.10e-3
G_AMPA['TuftedRS']['DeepLTS'] = 2.0e-3
G_NMDA['TuftedRS']['DeepLTS'] = 0.10e-3
G_AMPA['TuftedRS']['NontuftedRS'] = 1.0e-3
G_NMDA['TuftedRS']['NontuftedRS'] = 0.10e-3

#      G_GABA['DeepBasket']['SpinyStellate'] = 1.0e-3
G_GABA['DeepBasket']['SpinyStellate'] = 1.5e-3 # ? suppress spiny stellate bursts ?
G_GABA['DeepBasket']['TuftedIB'] = 0.7e-3
G_GABA['DeepBasket']['TuftedRS'] = 0.7e-3
G_GABA['DeepBasket']['DeepBasket'] = 0.2e-3
G_GABA['DeepBasket']['DeepAxoaxonic'] = 0.2e-3
G_GABA['DeepBasket']['DeepLTS'] = 0.7e-3
G_GABA['DeepBasket']['NontuftedRS'] = 0.7e-3

G_GABA['DeepAxoaxonic']['SupPyrRS'] = 1.0e-3
G_GABA['DeepAxoaxonic']['SupPyrFRB'] = 1.0e-3
#      G_GABA['DeepAxoaxonic']['SpinyStellate'] = 1.0e-3
G_GABA['DeepAxoaxonic']['SpinyStellate'] = 1.5e-3 # ? suppress spiny stellate bursts ?
G_GABA['DeepAxoaxonic']['TuftedIB'] = 1.0e-3
G_GABA['DeepAxoaxonic']['TuftedRS'] = 1.0e-3
G_GABA['DeepAxoaxonic']['NontuftedRS'] = 1.0e-3

G_GABA['DeepLTS']['SupPyrRS'] = .01e-3
G_GABA['DeepLTS']['SupPyrFRB'] = .01e-3
G_GABA['DeepLTS']['SupBasket'] = .01e-3
G_GABA['DeepLTS']['SupAxoaxonic'] = .01e-3
G_GABA['DeepLTS']['SupLTS'] = .05e-3
G_GABA['DeepLTS']['SpinyStellate'] = .01e-3
#      G_GABA['DeepLTS']['TuftedIB'] = .02e-3
G_GABA['DeepLTS']['TuftedIB'] = .05e-3 # will this help suppress bursting?
G_GABA['DeepLTS']['TuftedRS'] = .02e-3
G_GABA['DeepLTS']['DeepBasket'] = .01e-3
G_GABA['DeepLTS']['DeepAxoaxonic'] = .01e-3
G_GABA['DeepLTS']['DeepLTS'] = .05e-3
G_GABA['DeepLTS']['NontuftedRS'] = .01e-3

G_AMPA['TCR']['SupPyrRS'] = 0.5e-3
G_NMDA['TCR']['SupPyrRS'] = 0.05e-3
G_AMPA['TCR']['SupPyrFRB'] = 0.5e-3
G_NMDA['TCR']['SupPyrFRB'] = 0.05e-3
#      G_AMPA['TCR']['SupBasket'] = 1.0e-3
G_AMPA['TCR']['SupBasket'] = 0.1e-3
# try a variation in which main feedforward inhibtion from thalamus
# is via deep interneurons.  May be necessary later to include special
# layer 4 interneurons
#      G_NMDA['TCR']['SupBasket'] = .10e-3
G_NMDA['TCR']['SupBasket'] = .01e-3
#      G_AMPA['TCR']['SupAxoaxonic'] = 1.0e-3
G_AMPA['TCR']['SupAxoaxonic'] = 0.1e-3
#      G_NMDA['TCR']['SupAxoaxonic'] = .10e-3
G_NMDA['TCR']['SupAxoaxonic'] = .01e-3
G_AMPA['TCR']['SpinyStellate'] = 1.0e-3
G_NMDA['TCR']['SpinyStellate'] = .10e-3
G_AMPA['TCR']['TuftedIB'] = 1.5e-3
G_NMDA['TCR']['TuftedIB'] = .15e-3
G_AMPA['TCR']['TuftedRS'] = 1.5e-3
G_NMDA['TCR']['TuftedRS'] = .15e-3
#      G_AMPA['TCR']['DeepBasket'] = 1.0e-3
G_AMPA['TCR']['DeepBasket'] = 1.5e-3
G_NMDA['TCR']['DeepBasket'] = .10e-3
G_AMPA['TCR']['DeepAxoaxonic'] = 1.0e-3
G_NMDA['TCR']['DeepAxoaxonic'] = .10e-3
G_AMPA['TCR']['nRT'] = 0.75e-3   
G_NMDA['TCR']['nRT'] = .15e-3
G_AMPA['TCR']['NontuftedRS'] = 1.0e-3    
G_NMDA['TCR']['NontuftedRS'] = .10e-3

#      G_GABA['nRT']['TCR'] = 1.0e-3
# objref G_GABA['nRT']['TCR']
# 	G_GABA['nRT']['TCR'] = new Vector(num_nRT+1)

# The GABA conductance baseline is uniformly randomly distributed
# between 0.7 - 2.1 nS 
# (numpy.random.random_sample()*(2.1 - 0.7) + 0.7) * 1e-3 # nS -> uS
G_GABA['nRT']['TCR'] = 1.0 
# Values here need to be set below  
G_GABA['nRT']['nRT'] = 0.30e-3

G_AMPA['NontuftedRS']['SupPyrRS'] = 0.5e-3
G_NMDA['NontuftedRS']['SupPyrRS'] = 0.05e-3
G_AMPA['NontuftedRS']['SupPyrFRB'] = 0.5e-3
G_NMDA['NontuftedRS']['SupPyrFRB'] = 0.05e-3
G_AMPA['NontuftedRS']['SupBasket'] = 1.0e-3
G_NMDA['NontuftedRS']['SupBasket'] = 0.1e-3
G_AMPA['NontuftedRS']['SupAxoaxonic'] = 1.0e-3
G_NMDA['NontuftedRS']['SupAxoaxonic'] = 0.1e-3
G_AMPA['NontuftedRS']['SupLTS'] = 1.0e-3
G_NMDA['NontuftedRS']['SupLTS'] = 0.1e-3
G_AMPA['NontuftedRS']['SpinyStellate'] = 0.5e-3
G_NMDA['NontuftedRS']['SpinyStellate'] = 0.05e-3
G_AMPA['NontuftedRS']['TuftedIB'] = 1.0e-3
G_NMDA['NontuftedRS']['TuftedIB'] = 0.1e-3
G_AMPA['NontuftedRS']['TuftedRS'] = 1.0e-3
G_NMDA['NontuftedRS']['TuftedRS'] = 0.1e-3
G_AMPA['NontuftedRS']['DeepBasket'] = 3.0e-3
G_NMDA['NontuftedRS']['DeepBasket'] = .10e-3
G_AMPA['NontuftedRS']['DeepAxoaxonic'] = 3.0e-3
G_NMDA['NontuftedRS']['DeepAxoaxonic'] = .10e-3
G_AMPA['NontuftedRS']['DeepLTS'] = 2.0e-3
G_NMDA['NontuftedRS']['DeepLTS'] = .10e-3
G_AMPA['NontuftedRS']['TCR'] = .75e-3
G_NMDA['NontuftedRS']['TCR'] = .075e-3
G_AMPA['NontuftedRS']['nRT'] = 0.5e-3
G_NMDA['NontuftedRS']['nRT'] = 0.05e-3
G_AMPA['NontuftedRS']['NontuftedRS'] = 1.0e-3
G_NMDA['NontuftedRS']['NontuftedRS'] = 0.1e-3
# End defining synaptic conductance scaling factors

from config import ms, uS

# Here we update the units : from uS to S
for pre, gmap in G_AMPA.items():
    for post, g in gmap.items():
        G_AMPA[pre][post] = g * uS # uS to S
for pre, gmap in G_GABA.items():
    for post, g in gmap.items():
        G_GABA[pre][post] = g * uS # uS to S
for pre, gmap in G_NMDA.items():
    for post, g in gmap.items():
        G_NMDA[pre][post] = g * uS # uS to S

# Convert unit of time constant tau from ms to s
for pre, taumap in TAU_NMDA.items():
    for post, tau in taumap.items():
        TAU_NMDA[pre][post] = tau * ms # ms to s
for pre, taumap in TAU_AMPA.items():
    for post, tau in taumap.items():
        TAU_AMPA[pre][post] = tau * ms # ms to s
for pre, taumap in TAU_GABA.items():
    for post, tau in taumap.items():
        TAU_GABA[pre][post] = tau * ms # ms to s
for pre, taumap in TAU_GABA_FAST.items():
    for post, tau in taumap.items():
        TAU_GABA_FAST[pre][post] = tau * ms # ms to s
for pre, taumap in TAU_GABA_SLOW.items():
    for post, tau in taumap.items():
        TAU_GABA_SLOW[pre][post] = tau * ms # ms to s


SYNAPTIC_DELAY_DEFAULT = 0.05e-3
SYNAPTIC_DELAY_THALAMOCORTICAL = 1e-3
SYNAPTIC_DELAY_CORTICOTHALAMIC = 5e-3

PR = defaultdict(lambda : defaultdict(lambda: 1.0)) # Returns probability 1 for all undefined pairs.
PR['SpinyStellate']['SupPyrFRB'] = 0.79 # Silver RA, Lubke J, Sakmann B, Feldmeyer D (2003) Science 302:1981–1984.
PR['SpinyStellate']['SupPyrRS'] = 0.79 # Silver RA, Lubke J, Sakmann B, Feldmeyer D (2003) Science 302:1981–1984.

THRESHOLD_DEFAULT = 0.0



# 
# synapse.py ends here
