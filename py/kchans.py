# kchans.py --- 
# 
# Filename: kchans.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 17 23:58:49 2009 (+0530)
# Version: 
# Last-Updated: Wed Apr 22 00:06:36 2009 (+0530)
#           By: subhasis ray
#     Update #: 94
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

from channel import ChannelBase
from numpy import where, linspace, exp
import config

class KDR(ChannelBase):
    """Delayed rectifier current

    "In hippocampal pyramidal neurons, however, it has been reported have relatively slow activation, with a time to peak of some 50-100 msec and even slower inactivation. Such a slow activation would make it ill suited to participate in the repolarization of the AP.... An equation that can describe IK(DR) in cortical neurons is
    
    IK(DR) = m^3 * h * gbar_K(DR) * (Vm - EK)
    
    where m and h depend on voltage and time."
        - Johnston & Wu, Foundations of Cellular Neurophysiology (1995).

    But in Traub 2005, the equation used is:
    
    IK(DR) = m^4 * gbar_K(DR) * (Vm - EK)
    """

    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 4)
	self.Ek = -95e-3
	v = linspace(config.vmin, config.vmax, config.ndivs + 1) 
	tau_m = where(v < -10e-3, \
			   1e-3 * (0.25 + 4.35 * exp((v + 10.0e-3) / 10.0e-3)), \
			   1e-3 * (0.25 + 4.35 * exp((- v - 10.0e-3) / 10.0e-3)))
	m_inf = 1.0 / (1.0 + exp((- v - 29.5e-3) / 10e-3))
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
	self.xGate.tweakTau()

class KDR_FS(ChannelBase):
    """KDR for fast spiking neurons"""
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 4)
	self.Ek = -95e-3
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1.0 / (1.0 + exp((- v - 27e-3) / 11.5e-3))
	tau_m =  where(v < -10e-3, \
			   1e-3 * (0.25 + 4.35 * exp((v + 10.0e-3) / 10.0e-3)), \
			   1e-3 * (0.25 + 4.35 * exp((- v - 10.0e-3) / 10.0e-3)))
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = tau_m[i]
            self.xGate.B[i] = m_inf[i]
        self.xGate.tweakTau()
	
class KA(ChannelBase):
    """A type K+ channel"""
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 4, 1)
	self.Ek = -95e-3
        self.initX = 0.0
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

class K2(ChannelBase):
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 1, 1)
	self.Ek = -95e-3
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
	
class KM(ChannelBase):
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 1)
	self.Ek = -95e-3
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	a =  1e3 * 0.02 / ( 1 + exp((-v - 20e-3 ) / 5e-3))
	b = 1e3 * 0.01 * exp((-v - 43e-3) / 18e-3)
	for i in range(config.ndivs + 1):
            self.xGate.A[i] = a[i]
            self.xGate.B[i] = b[i]
	self.xGate.tweakAlpha()


	

# 
# kchans.py ends here
