# archan.py --- 
# 
# Filename: archan.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Mon Apr 27 15:34:07 2009 (+0530)
# Version: 
# Last-Updated: Sun May  3 23:30:18 2009 (+0530)
#           By: subhasis ray
#     Update #: 4
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

import config
import moose
from numpy import exp, linspace

from channel import ChannelBase

class AR(ChannelBase):
    """Combined cation current."""
    v = ChannelBase.v_array
    m_inf  = 1 / ( 1 + exp( ( v * 1e3 + 75 ) / 5.5 ) )
    tau_m = 1e-3 / ( exp( -14.6 - 0.086 * v * 1e3) + exp( -1.87 + 0.07 * v * 1e3))

    def __init__(self, name, parent, Ek=-35e-3):
	ChannelBase.__init__(self, name, parent, 1, 0)
	for i in range(len(self.xGate.A)):
	    self.xGate.A[i] = AR.tau_m[i]
	    self.xGate.B[i] = AR.m_inf[i]
	self.xGate.tweakTau()
	self.X = 0.25
	self.Ek = Ek



# 
# archan.py ends here
