# trb_globals.py --- 
# 
# Filename: trb_globals.py
# Description: global variables for PyMOOSE version of traub model
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Nov 29 03:27:46 2008 (+0530)
# Version: 
# Last-Updated: Tue Dec  9 09:21:23 2008 (+0530)
#           By: subhasis ray
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

# Code:

import moose

class Globals:
    """Contains global variables.

    These are essentially constant through out the simulation and
    should be set at the start of execution.
    """
    SIMULATOR = "pymoose"
    SIMDT = 1E-5
    PLOTDT = 1E-4
    VMIN = -120E-3
    VMAX = 40E-3
    NDIVS = 640
    dV = (VMAX - VMIN) / NDIVS
    INJECTION = 1E-11
    CONTEXT = moose.PyMooseBase.getContext()
    E_NA = 50e-3
    E_K = -95e-3
    E_K_FS = -100e-3
    E_CA = 125e-3
    E_AR = -35e-3
    # the following are set to be at the start of simulation
    simtime = 0.0
    plotsteps = 0
    simsteps = 0
    COMPARTMENT_FIELDS = ["len", "dia", "initVm", "Vm", "Em", "Rm", "Cm", "Ra", "inject", "Im"]
    HHCHANNEL_FIELDS = ["Ek", "Gbar", "Xpower", "Ypower", "Zpower", "instant", "Ik"]
# ! class Globals


# 
# trb_globals.py ends here
