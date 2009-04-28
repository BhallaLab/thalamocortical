# nachans.py --- 
# 
# Filename: nachans.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 17 23:58:13 2009 (+0530)
# Version: 
# Last-Updated: Wed Apr 29 02:19:19 2009 (+0530)
#           By: subhasis ray
#     Update #: 65
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

from numpy import where, linspace, exp
import pylab
import config
from channel import ChannelBase

class NaChannel(ChannelBase):
    """Dummy base class for all Na+ channels"""
    def __init__(self, name, parent, x, y=0):
        ChannelBase.__init__(self, name, parent, x, y)

class NaF(NaChannel):
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 3, 1)
        self.Ek = 50e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1)
        tau_m = where(v < -30e-3, \
                          1.0e-3 * (0.025 + 0.14 * exp((v + 30.0e-3) / 10.0e-3)), \
                          1.0e-3 * (0.02 + 0.145 * exp(( - v - 30.0e-3) / 10.0e-3)))

        tau_h = 1.0e-3 * (0.15 + 1.15 / ( 1.0 + exp(( v + 37.0e-3) / 15.0e-3)))
        
        m_inf = 1.0 / (1.0 + exp(( - v - 38e-3) / 10e-3))
        
        h_inf = 1.0 / (1.0 + exp((v + 62.9e-3) / 10.7e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
        self.yGate.tweakTau()
        self.xGate.A.dumpFile("naf_xa.plot")
        self.xGate.B.dumpFile("naf_xb.plot")
        self.yGate.A.dumpFile("naf_ya.plot")
        self.yGate.B.dumpFile("naf_yb.plot")
        
class NaF2(NaChannel):
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 3, 1)
        self.Ek = 50e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1)
        tau_m = where(v < -30e-3, \
                          1.0e-3 * (0.0125 + 0.1525 * exp ((v + 30e-3) / 10e-3)), \
                          1.0e-3 * (0.02 + 0.145 * exp((-v - 30e-3) / 10e-3)))
        
        tau_h = 1e-3 * (0.225 + 1.125 / ( 1 + exp( (  v  + 37e-3 ) / 15e-3 ) ))
        
        h_inf = 1.0 / (1.0 + exp((v + 58.3e-3) / 6.7e-3))

        m_inf = 1.0 / (1.0 + exp(( - v - 38e-3) / 10e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
        self.yGate.tweakTau()

class NaP(NaChannel):
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 1)
        self.Ek = 50e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1)
        tau_m = where(v < -40e-3, \
                          1.0e-3 * (0.025 + 0.14 * exp((v + 40e-3) / 10e-3)), \
                          1.0e-3 * (0.02 + 0.145 * exp((-v - 40e-3) / 10e-3)))
        m_inf = 1.0 / (1.0 + exp((-v - 48e-3) / 10e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()
        self.xGate.A.dumpFile("nap_xa.plot")
        self.xGate.B.dumpFile("nap_xb.plot")


class NaPF(NaChannel):
    """Persistent Na+ current, fast"""
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 3)
        self.Ek = 50e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1)
        tau_m = where(v < -30e-3, \
                           1.0e-3 * (0.025 + 0.14 * exp((v  + 30.0e-3) / 10.0e-3)), \
                           1.0e-3 * (0.02 + 0.145 * exp((- v - 30.0e-3) / 10.0e-3)))
        m_inf = 1.0 / (1.0 + exp((-v - 38e-3) / 10e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()

class NaPF_TCR(NaChannel):
    """Persistent Na+ channel specific to TCR cells. Only difference
    with NaPF is power of m is 1 as opposed 3."""
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 1)
        self.Ek = 50e-3
        shift = 7e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1) + shift
        tau_m = where(v < -30e-3, \
                           1.0e-3 * (0.025 + 0.14 * exp((v  + 30.0e-3) / 10.0e-3)), \
                           1.0e-3 * (0.02 + 0.145 * exp((- v - 30.0e-3) / 10.0e-3)))
        m_inf = 1.0 / (1.0 + exp((-v - 38e-3) / 10e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()

class NaF_TCR(NaChannel):
    """Fast Na+ channel for TCR cells. This is almost identical to
    NaF, but there is a nasty voltage shift in the tables."""
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 3, 1)
        self.Ek = 50e-3
        shift_mnaf = -5.5e-3
        shift_hnaf = -7e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1) 
        tau_h = 1.0e-3 * (0.15 + 1.15 / ( 1.0 + exp(( v + 37.0e-3) / 15.0e-3)))        
        v = v + shift_hnaf
        h_inf = 1.0 / (1.0 + exp((v + 62.9e-3) / 10.7e-3))
        v = v - shift_hnaf + shift_mnaf
        tau_m = where(v < -30e-3, \
                          1.0e-3 * (0.025 + 0.14 * exp((v + 30.0e-3) / 10.0e-3)), \
                          1.0e-3 * (0.02 + 0.145 * exp(( - v - 30.0e-3) / 10.0e-3)))
        m_inf = 1.0 / (1.0 + exp(( - v - 38e-3) / 10e-3))

        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
            self.yGate.A[i] = tau_h[i]
            self.yGate.B[i] = h_inf[i]
        self.xGate.tweakTau()
        self.yGate.tweakTau()
        self.xGate.A.dumpFile("naf_xa.plot")
        self.xGate.B.dumpFile("naf_xb.plot")
        self.yGate.A.dumpFile("naf_ya.plot")
        self.yGate.B.dumpFile("naf_yb.plot")
        
class NaPF_SS(NaChannel):
    def __init__(self, name, parent):
        NaChannel.__init__(self, name, parent, 3)
        self.Ek = 50e-3
        shift = -2.5e-3
        v = linspace(config.vmin, config.vmax, config.ndivs + 1) + shift
        tau_m = where(v < -30e-3, \
                           1.0e-3 * (0.025 + 0.14 * exp((v  + 30.0e-3) / 10.0e-3)), \
                           1.0e-3 * (0.02 + 0.145 * exp((- v - 30.0e-3) / 10.0e-3)))
        m_inf = 1.0 / (1.0 + exp((- v - 38e-3) / 10e-3))
        for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()

# 
# nachans.py ends here
