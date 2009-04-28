# kchans.py --- 
# 
# Filename: kchans.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 17 23:58:49 2009 (+0530)
# Version: 
# Last-Updated: Wed Apr 29 02:18:58 2009 (+0530)
#           By: subhasis ray
#     Update #: 391
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
from channel import ChannelBase
from numpy import where, linspace, exp, arange, ones, zeros
import config


class KChannel(ChannelBase):
    """This is a dummy base class to keep type information."""
    def __init__(self, name, parent, xpower=1, ypower=0):
        ChannelBase.__init__(self, name, parent, xpower, ypower)
        self.Ek = -95e-3


class KDR(KChannel):
    """Delayed rectifier current

    "In hippocampal pyramidal neurons, however, it has been reported have relatively slow activation, with a time to peak of some 50-100 msec and even slower inactivation. Such a slow activation would make it ill suited to participate in the repolarization of the AP.... An equation that can describe IK(DR) in cortical neurons is
    
    IK(DR) = m^3 * h * gbar_K(DR) * (Vm - EK)
    
    where m and h depend on voltage and time."
        - Johnston & Wu, Foundations of Cellular Neurophysiology (1995).

    But in Traub 2005, the equation used is:
    
    IK(DR) = m^4 * gbar_K(DR) * (Vm - EK)
    """

    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 4)
	v = linspace(config.vmin, config.vmax, config.ndivs + 1) 
	tau_m = where(v < -10e-3, \
			   1e-3 * (0.25 + 4.35 * exp((v + 10.0e-3) / 10.0e-3)), \
			   1e-3 * (0.25 + 4.35 * exp((- v - 10.0e-3) / 10.0e-3)))
	m_inf = 1.0 / (1.0 + exp((- v - 29.5e-3) / 10e-3))
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
	self.xGate.tweakTau()


class KDR_FS(KChannel):
    """KDR for fast spiking neurons"""
    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 4)
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1.0 / (1.0 + exp((- v - 27e-3) / 11.5e-3))
	tau_m =  where(v < -10e-3, \
			   1e-3 * (0.25 + 4.35 * exp((v + 10.0e-3) / 10.0e-3)), \
			   1e-3 * (0.25 + 4.35 * exp((- v - 10.0e-3) / 10.0e-3)))
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()
	

class KA(KChannel):
    """A type K+ channel"""
    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 4, 1)
        self.X = 0.0
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1 / ( 1 + exp( ( - v - 60e-3 ) / 8.5e-3 ) )
	tau_m =  1e-3 * (0.185 + 0.5 / ( exp( ( v + 35.8e-3 ) / 19.7e-3 ) + exp( ( - v - 79.7e-3 ) / 12.7e-3 ) ))
	h_inf =   1 / ( 1 + exp( ( v + 78e-3 ) / 6e-3 ) )
	tau_h = where( v <= -63e-3,\
                           1e-3 * 0.5 / ( exp( ( v + 46e-3 ) / 5e-3 ) + exp( ( - v - 238e-3 ) / 37.5e-3 ) ), \
                           9.5e-3)

	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
	self.yGate.tweakTau()


class KA_IB(KChannel):
    """A type K+ channel for tufted intrinsically bursting cells -
    multiplies tau_h of KA by 2.6"""
    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 4, 1)
        self.X = 0.0
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1 / ( 1 + exp( ( - v - 60e-3 ) / 8.5e-3 ) )
	tau_m =  1e-3 * (0.185 + 0.5 / ( exp( ( v + 35.8e-3 ) / 19.7e-3 ) + exp( ( - v - 79.7e-3 ) / 12.7e-3 ) ))
	h_inf =   1 / ( 1 + exp( ( v + 78e-3 ) / 6e-3 ) )
	tau_h = 2.6 * where( v <= -63e-3,\
                                 1e-3 * 0.5 / ( exp( ( v + 46e-3 ) / 5e-3 ) + exp( ( - v - 238e-3 ) / 37.5e-3 ) ), \
                                 9.5e-3)

	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
	self.yGate.tweakTau()


class K2(KChannel):
    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 1, 1)
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1.0 / (1 + exp((-v - 10e-3) / 17e-3))
	tau_m = 1e-3 * (4.95 + 0.5 / (exp((v - 81e-3) / 25.6e-3) + \
					  exp((-v - 132e-3) / 18e-3)))
	
	h_inf = 1.0 / (1 + exp((v + 58e-3) / 10.6e-3))
	tau_h = 1e-3 * (60 + 0.5 / (exp((v - 1.33e-3) / 200e-3) + \
					exp((-v - 130e-3) / 7.1e-3)))
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
	self.yGate.tweakTau()
	

class KM(KChannel):
    def __init__(self, name, parent):
	KChannel.__init__(self, name, parent, 1)
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	a =  1e3 * 0.02 / ( 1 + exp((-v - 20e-3 ) / 5e-3))
	b = 1e3 * 0.01 * exp((-v - 43e-3) / 18e-3)
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = a[i]
            self.xGate.B[i] = b[i]
	self.xGate.tweakAlpha()


class KCaChannel(KChannel):
    """[Ca+2] dependent K+ channel base class."""
    def __init__(self, name, parent, xpower=0, ypower=0, zpower=1):
        KChannel.__init__(self, name, parent, xpower, ypower)
        self.connected_to_ca = False
        self.Zpower = zpower
        self.zGate = moose.HHGate('zGate', self)
        self.zGate.A.xmin = 0.0
        self.zGate.A.xdivs = 1000
        self.zGate.B.xmin = 0.0
        self.zGate.B.xdivs = 1000        
        self.zGate.A.calcMode = 1
        self.zGate.B.calcMode = 1
        self.useConcentration = True


class KAHPBase(KCaChannel):
    def __init__(self, name, parent):
        KCaChannel.__init__(self, name, parent)
        self.zGate.A.xmax = 1.0
        self.zGate.B.xmax = 1.0
        

class KAHP(KAHPBase):
    """AHP type K+ current"""
    def __init__(self, name, parent):
        KAHPBase.__init__(self, name, parent)

        ca_conc = linspace(self.zGate.A.xmin, self.zGate.A.xmax, self.zGate.A.xdivs + 1)
        alpha = where(ca_conc < 100.0 * 1e-3 , 0.1 * ca_conc, 10.0)
        beta =  ones(self.zGate.B.xdivs + 1) * 10.0
        for i in range(len(alpha)):
            self.zGate.A[i] = alpha[i]
            self.zGate.B[i] = beta[i]
        self.zGate.tweakAlpha()



class KAHP_SLOWER(KAHPBase):
    def __init__(self, name, parent):
        KAHPBase.__init__(self, name, parent)
        ca_conc = linspace(self.zGate.A.xmin, self.zGate.A.xmax, self.zGate.A.xdivs + 1)
        alpha = where(ca_conc < 500.0e-3, 1e6 * ca_conc / 50000, 10.0)
        beta =  ones(self.zGate.B.xdivs + 1) * 1.0
        for i in range(len(alpha)):
            self.zGate.A[i] = alpha[i]
            self.zGate.B[i] = beta[i]
        self.zGate.tweakAlpha()


class KAHP_DP(KAHPBase):
    """KAHP for deep pyramidal cell"""
    def __init__(self, name, parent):
        KAHPBase.__init__(self, name, parent)
        ca_conc = linspace(self.zGate.A.xmin, self.zGate.A.xmax, self.zGate.A.xdivs + 1)
        alpha = where(ca_conc < 100.0 * 1e-3, 1e-4 * ca_conc, 0.01)
        beta =  0.001
        for i in range(len(alpha)):
            self.zGate.A[i] = alpha[i]
            self.zGate.B[i] = beta[i]
        self.zGate.tweakAlpha()

class KC(KCaChannel):
    """C type K+ channel"""
    def __init__(self, name, parent):
        KCaChannel.__init__(self, name, parent, 1, 0, 1)
        self.zGate.A.xmax = 1.0
        self.zGate.B.xmax = 1.0
        ca_conc = linspace(self.zGate.A.xmin, self.zGate.A.xmax, self.zGate.A.xdivs + 1)
        alpha = where(ca_conc < 250e-3, ca_conc / 250e-3, 1.0)
        for i in range(len(alpha)):
            self.zGate.A[i] = alpha[i]
            self.zGate.B[i] = 1.0
        self.instant = 4 # Zgate m is instantaneous: m = A/B
        self.zGate.A.calcMode = 1
        self.zGate.B.calcMode = 1
        self.X = 0.0
        v = linspace(config.vmin, config.vmax, config.ndivs + 1)
        alpha = zeros(len(v))
        beta = zeros(len(v))
        alpha = where(v < -10e-3, 
                      2e3 / 37.95 * ( exp( ( v * 1e3 + 50 ) / 11 - ( v * 1e3 + 53.5 ) / 27 ) ),
                      2e3 * exp(( - v * 1e3- 53.5) / 27))
        beta = where( v < -10e-3,
                      2e3 * exp(( - v * 1e3 - 53.5) / 27) - alpha, 
                      0.0)
        for i in range(len(alpha)):
            self.xGate.A[i] = alpha[i]
            self.xGate.B[i] = beta[i]
        self.xGate.tweakAlpha()
        self.X = 0.0

        
class KC_FAST(KC):
    """Fast KC channel"""
    def __init__(self, name, parent):
        KC.__init__(self, name, parent)
        for i in range(len(self.xGate.A)):
            self.xGate.A[i] = 2 * self.xGate.A[i]
            self.xGate.B[i] = 2 * self.xGate.B[i]
        
if __name__ == "__main__":
    a = KC_FAST('kc', moose.Neutral('/'))
# 
# kchans.py ends here
