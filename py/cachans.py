# cachans.py --- 
# 
# Filename: cachans.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Sat Apr 18 00:18:24 2009 (+0530)
# Version: 
# Last-Updated: Sat Apr 18 01:04:59 2009 (+0530)
#           By: subhasis ray
#     Update #: 46
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

class CaL(ChannelBase):
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 2)
	self.Ek = 125e-3
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	alpha = 1.6e3 / (1.0 + exp(-0.072e3 * (v - 5e-3)))
	beta = where( v < (1e-6 - 8.9e-3),
		      100 * exp(-v / 5),
		      20 * v / (exp(v / 5) - 1))
	for i in range(config.ndivs + 1):
	    self.xGate.A[i] = alpha[i]
	    self.xGate.B[i] = beta[i]
	self.xGate.tweakAlpha()

class CaT(ChannelBase):
    def __init__(self, name, parent):
	ChannelBase.__init__(self, name, parent, 2, 1)
	self.Ek = 125e-3
	v = linspace(config.vmin, config.vmax, config.ndivs + 1)
	m_inf = 1 / (1 + exp( (- v - 56e-3) / 6.2e-3))
	tau_m = 1e-3 * (0.204 + 0.333 / ( exp(( v + 15.8e-3) / 18.2e-3 ) + 
					  exp((- v - 131e-3) / 16.7e-3)))
	h_inf = 1 / (1 + exp(( v + 80e-3 ) / 4e-3))
	tau_h = where( v < -81e-3, 
		       1e-3 * 0.333 * exp( ( v + 466e-3 ) / 66.6e-3 ),
		       1e-3 * (9.32 + 0.333 * exp( ( -v - 21 ) / 10.5 )))
	for i in range(config.ndivs + 1):
	    self.xGate.A[i] = tau_m[i]
	    self.xGate.B[i] = m_inf[i]
	    self.yGate.A[i] = tau_h[i]
	    self.yGate.B[i] = h_inf[i]
	self.xGate.tweakTau()
	self.yGate.tweakTau()
		       



# 
# cachans.py ends here
