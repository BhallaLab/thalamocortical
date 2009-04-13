# trb_channels.py --- 
# 
# Filename: trb_channels.py
# Description: Implemetation of Na channels.
# Author: Subhasis Ray
# Maintainer: 
# Created: Wed Dec  3 13:56:20 2008 (+0530)
# Version: 
# Last-Updated: Sat Dec 20 11:05:12 2008 (+0530)
#           By: subhasis ray
#     Update #: 349
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
# 2008-12-03: combined old code from naf.py and nap.py.
# 
# 
# 


# Code:

from math import *

import moose

from trb_globals import *
import trb_utility as util

class NaF(moose.HHChannel):
    """Traub 2005: fast Na channels"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, yGate, = util.init_changates(self, 3, 1)
        self.Ek = Globals.E_NA
        v = Globals.VMIN
        for i in range(0,Globals.NDIVS+1):
            xGate.A[i] = self.calc_tau_m(v)
            xGate.B[i] = self.calc_m_inf(v)
            yGate.A[i] = self.calc_tau_h(v)
            yGate.B[i] = self.calc_h_inf(v)
            v += Globals.dV
        xGate.tweakTau()
        yGate.tweakTau()
#         if not Globals.CONTEXT.connect(self.parent, "channel", self.id, "channel"):
#             logging.error("error connecting %s to %s", self.parent.path, self.path)

    def calc_m_inf(self, v):
        return 1.0/(1.0+ exp((-v-38e-3)/10e-3))

    def calc_h_inf(self, v):
        return 1.0/(1.0+ exp((v+62.9e-3)/10.7e-3))

    def calc_tau_m(self, v):        
        if v < -30e-3:
            return 1.0e-3*(0.025 + 0.14 * exp((v+30.0e-3)/10.0e-3))
        else:
            return 1.0e-3*(0.02 + 0.145 * exp((-v-30.0e-3)/10.0e-3))

    def calc_tau_h(self, v):
        return 1.0e-3*(0.15 + 1.15 / ( 1.0 + exp((v+37.0e-3)/15.0e-3)))
      
class NaF2(NaF):
    """Fast Na channel present in most of the cell types in Traub et al 2005"""
    def __init__(self, *args):
        NaF.__init__(self, *args)
        xGate = moose.HHGate(self.path + "/xGate")
        yGate = moose.HHGate(self.path + "/yGate")
        xGate.A.dumpFile("pymoose_xa_" + self.name + ".plot")
        xGate.B.dumpFile("pymoose_xb_" + self.name + ".plot")
        yGate.A.dumpFile("pymoose_ya_" + self.name + ".plot")
        yGate.A.dumpFile("pymoose_yb_" + self.name + ".plot")
    
    def calc_h_inf(self, v):
        return 1.0 / (1.0 + exp((v + 58.3e-3) / 6.7e-3))

    def calc_tau_m(self, v):
        if v < -30e-3:
            return 1.0e-3 * (0.0125 + 0.1525 * exp ((v + 30e-3) / 10e-3))
        else:
            return 1.0e-3 * (0.02 + 0.145 * exp((-v - 30e-3) / 10e-3))

    def calc_tau_h(self, v):
        return 1e-3 * (0.225 + 1.125 / ( 1 + exp( (  v  + 37e-3 ) / 15e-3 ) ))
        
# ! NaF2


class NaP(moose.HHChannel):
    """Persistent Na channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, = util.init_changates(self, 1, 0)
        self.Ek = Globals.E_NA
        v = Globals.VMIN
        for i in range(0,Globals.NDIVS+1):
            xGate.A[i] = self.calc_tau_m(v)
            xGate.B[i] = self.calc_m_inf(v)
            v += Globals.dV
        xGate.tweakTau()

    def calc_m_inf(self, v):
        return 1.0/(1.0+ exp((-v-40e-3)/10e-3))

    def calc_tau_m(self, v):        
        if v < -40e-3:
            return 1.0e-3 * (0.025 + 0.14 * exp((v + 40e-3) / 10e-3))
        else:
            return 1.0e-3 * (0.02 + 0.145 * exp((-v - 40e-3) / 10e-3))

class NaPTCR(NaP):
    """Persistent Na channel for TCR cells. This is identical to
NaFastGluta with a shift of 7 mv and conductance is linearly related
to activation variable"""
    def __init__(self, *args):
        apply(NaP.__init__, (self,)+args)

    def calc_m_inf(self, v):
        return 1.0/(1.0+ exp((-v-7e-3- 38e-3)/10e-3))

    def calc_tau_m(self, v):        
        if (v+7e-3) < -30e-3:
            return 1.0e-3*(0.025 + 0.14 * exp((v + 7e-3 + 30.0e-3)/10.0e-3))
        else:
            return 1.0e-3*(0.02 + 0.145 * exp((-v - 7e-3 - 30.0e-3)/10.0e-3))
    
class NaPSS(NaP):
    """Persistent Na+ current. Identical to NaFast except that here there
is no inactivation variable and the resting potential is shifted by -2.5 mV."""
    def __init__(self, *args):
        NaP.__init__(self, *args)
        self.Xpower = 3.0

    def calc_m_inf(self, v):
        return 1.0/(1.0+ exp((-v + 2.5e-3 - 38e-3)/10e-3))

    def calc_tau_m(self, v):        
        if (v-2.5e-3) < -30e-3:
            return 1.0e-3*(0.025 + 0.14 * exp((v - 2.5e-3 + 30.0e-3)/10.0e-3))
        else:
            return 1.0e-3*(0.02 + 0.145 * exp((-v + 2.5e-3 - 30.0e-3)/10.0e-3))


class NaPF(NaP):
    """Persistent Na+ current, fast"""
    def __init__(self, *args):
        apply(NaP.__init__, (self,)+args)
        self.Xpower = 3.0

    def calc_m_inf(self, v):
        return 1.0/(1.0+ exp((-v - 38e-3)/10e-3))

    def calc_tau_m(self, v):        
        if v < -30e-3:
            return 1.0e-3*(0.025 + 0.14 * exp((v  + 30.0e-3)/10.0e-3))
        else:
            return 1.0e-3*(0.02 + 0.145 * exp((- v - 30.0e-3)/10.0e-3))


class KDR(moose.HHChannel):
    """Delayed rectifier current

    "In hippocampal pyramidal neurons, however, it has been reported have relatively slow activation, with a time to peak of some 50-100 msec and even slower inactivation. Such a slow activation would make it ill suited to participate in the repolarization of the AP.... An equation that can describe IK(DR) in cortical neurons is
    
    IK(DR) = m^3 * h * gbar_K(DR) * (Vm - EK)
    
    where m and h depend on voltage and time."
        - Johnston & Wu, Foundations of Cellular Neurophysiology (1995).

    But in Traub 2005, the equation used is:
    
    IK(DR) = m^4 * gbar_K(DR) * (Vm - EK)
    """
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        self.Ek = Globals.E_K
        xGate, = util.init_changates(self, 4)
        v = Globals.VMIN
        for ii in range(0, Globals.NDIVS + 1):
            xGate.A[ii] = self.calc_tau_m(v)
            xGate.B[ii] = self.calc_m_inf(v)
            v += Globals.dV
        xGate.tweakTau()
#         if not Globals.CONTEXT.connect(self.parent, "channel", self.id, "channel"):
#             logging.error("error connecting %s to %s", self.parent.path, self.path)

    def calc_m_inf(self, v):
        return 1.0 / (1.0 + exp((- v - 29.5e-3) / 10e-3))

    def calc_tau_m(self, v):        
        if v < -10e-3:
            return 1e-3 * (0.25 + 4.35 * exp((v + 10.0e-3) / 10.0e-3))
        else:
            return 1e-3 * (0.25 + 4.35 * exp((- v - 10.0e-3) / 10.0e-3))
    # !calc_tau_m

#! class KDR


class KDRFS(KDR):
    """KDR for fast spiking neurons."""
    def __init__(self, *args):
        KDR.__init__(self, *args)
        self.Ek = Globals.E_K_FS

    def calc_m_inf(self, v):
        """Overrides the function with same name in KDR"""
        return 1.0 / (1.0 + exp((- v - 27e-3) / 11.5e-3))
#! class KDRFS


class KA(moose.HHChannel):
    """A type K+ channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, yGate, = util.init_changates(self, 4, 1)
        v = Globals.VMIN
        for i in range(0, Globals.NDIVS+1):
            xGate.A[i] = self.calc_tau_m(v)
            xGate.B[i] = self.calc_m_inf(v)
            yGate.A[i] = self.calc_tau_h(v)
            yGate.B[i] = self.calc_h_inf(v)
            v += Globals.dV
        xGate.tweakTau()
        yGate.tweakTau()
    # ! __init__ 

    def calc_m_inf(self, v):
        return 1.0 / (1.0 + exp((-v - 60e-3) / 8.5e-3))

    def calc_tau_m(self, v):        
        return 1.0e-3 * (0.185 + 0.5 /( exp ((v + 35.8e-3) / 19.7e-3) + exp((-v - 79.7e-3) / 12.7e-3)))

    def calc_tau_h( self, v ):
        if ( v <= -63e-3):
            return  1e-3 * 0.5/ ( exp ( (v + 46e-3) / 5e-3 ) + exp ( ( -v - 238e-3 )/37.5e-3 ) ) 
        else:
            return 9.5e-3

    def calc_h_inf( self, v ):
        return  1.0 / (1.0 + exp ( (v + 78e-3) / 6e-3 ) )


# ! class KA        


class K2(moose.HHChannel):
    """K2 channel"""
    def __init__(self, *args):
        apply(moose.HHChannel.__init__, (self,) + args)
        (xGate, yGate,) = util.init_changates(self, 1, 1)
        v = Globals.VMIN
        for i in range(0, (Globals.NDIVS + 1)):
            xGate.A[i] = self.calc_tau_m(v)
            xGate.B[i] = self.calc_m_inf(v)
            yGate.A[i] = self.calc_tau_h(v)
            yGate.B[i] = self.calc_h_inf(v)
            v += Globals.dV

        xGate.tweakTau()
        yGate.tweakTau()

    def calc_m_inf(self, v):
        return 1.0 / (1 + exp((-v - 10e-3) / 17e-3))

    def calc_tau_m(self, v):
        return 1e-3 * (4.95 + 0.5 / (exp((v - 81e-3) / 25.6e-3) + exp((-v - 132e-3) / 18e-3)))

    def calc_tau_h(self, v):
        return 1e-3 * (60 + 0.5 / (exp((v - 1.33e-3) / 200e-3) + exp((-v - 130e-3) / 7.1e-3)))

    def calc_h_inf(self, v):
        return 1.0 / (1 + exp((v + 58e-3) / 10.6e-3))
# ! class K2

class KM(moose.HHChannel):
    """Mascarinic sensitive K+ channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, = util.init_changates(self, 1, 0)
        v = Globals.VMIN
        for i in range(0, Globals.NDIVS + 1):
            xGate.A[i] = 0.02 / ( 1 + exp((-v - 20e-3 ) / 5e-3))
            xGate.B[i] = 0.01 * exp((-v - 43e-3) / 18e-3)
            v = v + Globals.dV
        xGate.tweakAlpha()
# !class KM

class AR(moose.HHChannel):
    """AR channel - combined current from miscellaneous cations"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, = util.init_changates(self, 1, 0)
        v = Globals.VMIN
        for i in range(0, Globals.NDIVS + 1):
            xGate.A[i] = 1e-3 / exp((v + 75e-3) / 5.5e-3)
            xGate.B[i] = 1e-3 / (exp( -14.6 - 86 * v ) + exp( - 1.87 + 70 * v))
            v = v + Globals.dV
        xGate.tweakTau()
# ! class AR

class KC(moose.HHChannel):
    """KC channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        self.useConcentration = True
        xGate, = util.init_changates(self, 1, 0)
        v = Globals.VMIN
        alpha = 0.0
        beta = 0.0
        for i in range(0, Globals.NDIVS + 1):
            if v < -10e-3:
                alpha =  (2e3 / 37.95) * exp ((v + 50e-3) / 11e-3 - (v + 53.5e-3) / 27e-3)
                beta = 2e3 * exp( ( - v - 53.5e-3 ) / 27e-3 ) - alpha
            else:
                alpha = 2e3 * exp( (- v - 53.5e-3) / 27e-3)
                beta = 0.0
            xGate.A[i] = alpha
            xGate.B[i] = beta
            v = v + Globals.dV
        xGate.tweakAlpha()
        xGate.A.calcMode = 0
        xGate.B.calcMode = 0
        # setup the z gate ( [Ca+2] dependency )
        xmin = 0.0 
        xmax = 250.0
        xdivs = 3000
        self.Zpower = 1
        zGate = moose.HHGate(self.path + "/zGate")
        zGate.A.xmin = xmin
        zGate.B.xmin = xmin
        zGate.A.xmax = xmax
        zGate.B.xmax = xmax
        zGate.A.xdivs = xdivs
        zGate.B.xdivs = xdivs
        
        chi = xmin
        dchi = (xmax - xmin) / xdivs
        for i in range(0, xdivs):
            zGate.A[i] = chi / 250.0
            zGate.B[i] = 1.0
            chi = chi + dchi
        # for any index above the range, Interpol returns the last value
        zGate.A[xdivs] = 1.0
        zGate.B[xdivs] = 1.0
        zGate.A.calcMode = 0
        zGate.B.calcMode = 0
        self.instant = 4 # INSTANTZ is 4
#         zGate.tabFill(3000, 2)
#!class KC


class KCF(KC):
    """KC channel with faster kinetics"""
    def __init__(self, *args):
        KC.__init__(self, *args)
        xGate = moose.HHGate(self.path + "/xGate")
        xdivs = self.xGate.A.xdivs
        for i in range(0, xdivs):
            xGate.A[i] = xGate.A[i] * 2
            xGate.B[i] = xGate.B[i] * 2
#! class KCF


class KAHP(moose.HHChannel):
    """After-hyperpolarization current"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        self.useConcentration = True
        xmin = -10.0 # arbitrary
        xmax = 1e6 # as large as possible
        xdivs = 100 # only upto [chi]
        self.Zpower = 1
        zGate = moose.HHGate(self.path + "/zGate")
        zGate.A.xmin = xmin
        zGate.B.xmin = xmin
        zGate.A.xmax = xmax
        zGate.B.xmax = xmax
        zGate.A.xdivs = xdivs
        zGate.B.xdivs = xdivs
        
        chi = xmin
        dchi = 1.0
        for i in range(0, xdivs):
            zGate.A[i] = chi / 10.0
            zGate.B[i] = 10.0
            chi = chi + dchi
        # for any index above the range, Interpol returns the last value
        zGate.A[xdivs] = 10.0
        zGate.B[xdivs] = 10.0
        zGate.A.calcMode = 0
        zGate.B.calcMode = 0
        zGate.tweakAlpha()
#! class KAHP

class KAHPSlow(moose.HHChannel):
    """Slow KAHP channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        self.useConcentration = True
        xmin = 0.0
        xmax = 500.0
        xdivs = 3000 
        self.Zpower = 1
        zGate = moose.HHGate(self.path + "/zGate")
        zGate.A.xmin = xmin
        zGate.B.xmin = xmin
        zGate.A.xmax = xmax
        zGate.B.xmax = xmax
        zGate.A.xdivs = xdivs
        zGate.B.xdivs = xdivs
        
        chi = xmin
        dchi = (xmax - xmin) / xdivs
        for i in range(0, xdivs):
            zGate.A[i] = chi * 0.02
            zGate.B[i] = 10.0
            chi = chi + dchi
        # for any index above the range, Interpol returns the last value
        zGate.A[xdivs] = 10.0
        zGate.B[xdivs] = 10.0
        zGate.A.calcMode = 0
        zGate.B.calcMode = 0
        zGate.tweakAlpha()
#! class KAHPSlow        


class CaL(moose.HHChannel):
    """Long lasting Ca2+ current"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        xGate, = util.init_changates(self, 2.0, 0.0)
        v = Globals.VMIN + 8.9e-3
        for i in range(0, Globals.NDIVS + 1):
            xGate.A[i] = 1.6e3 / ( 1.0 + exp( -0.072e3 * (v - 8.9e-3 - 5e-3)))
            if abs(v) < 1e-9:
                xGate.B[i] = 100 * exp (-v / 5)
            else:
                xGate.B[i] = 20 * v / (exp(v / 5) - 1)
            v = v + Globals.dV
        xGate.tweakAlpha()
#! class CaL


class CaT(moose.HHChannel):
    """Transient low-threshold Ca2+ current"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        (xGate, yGate) = util.init_changates(self, 2.0, 1.0)
        v = Globals.VMIN
        for i in range(0, Globals.NDIVS + 1):
            xGate.A[i] = 1e-3 * ( 0.204 + 0.333 / (exp((v + 15.8e-3) / 18.2e-3) + exp((-v - 131e-3) / 16.7e-3)))
            xGate.B[i] = 1.0 / (1 + exp((-v - 56e-3) / 6.2e-3))
            if v < -81e-3:
                yGate.A[i] = 1e-3 * 0.333 * exp((v + 466e-3) / 66.6e-3)
            else:
                yGate.A[i] = 1e-3 * (9.32 + 0.333 * exp((-v - 21e-3) / 10.5e-3))
            yGate.B[i] = 1.0 / ( 1 + exp(( v + 80e-3) / 4e-3))
            v = v + Globals.dV
        xGate.tweakTau()
        yGate.tweakTau()
#! class CaT


class CaTA(moose.HHChannel):
    """A variant of transient Ca2+ channel"""
    def __init__(self, *args):
        moose.HHChannel.__init__(self, *args)
        (xGate, yGate) = util.init_changates(self, 2.0, 1.0)
        v = Globals.VMIN
        for i in range(0, Globals.NDIVS + 1):
            xGate.A[i] = 1e-3 * (1 + 0.33 / (exp ((v + 27e-3) / 10e-3) + exp ((-v - 102e-3) / 15e-3)))
            xGate.B[i] = 1.0 / (1.0 + exp ((-v - 52e-3) / 7.4e-3))
            yGate.A[i] = 1e-3 * (28.3 + 0.33 / (exp ((v + 48e-3) / 4e-3) + exp ((-v - 407e-3) / 50e-3)))
            yGate.B[i] = 1.0 / (1 + exp ((v + 80e-3) / 5e-3))
            v = v + Globals.dV
        xGate.tweakTau()
        yGate.tweakTau()
# !class CaTA


# 
# trb_channels.py ends here
