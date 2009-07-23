# config.py --- 
# 
# Filename: config.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 17 14:36:30 2009 (+0530)
# Version: 
# Last-Updated: Thu Jul 23 22:40:57 2009 (+0530)
#           By: subhasis ray
#     Update #: 22
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This contains global variables used in other modules.
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

import sys
sys.path.append('/home/src/sim/cortical/py')
import moose
context = moose.PyMooseBase.getContext()
lib = moose.Neutral('/library')
root = moose.Neutral("/")

simdt = 1e-5
plotdt = 1e-5
vmin = -120e-3
vmax = 40e-3
ndivs = 640
dv = (vmax - vmin) / ndivs
channel_name_list = ['AR','CaPool','CaL','CaT','CaT_A','K2','KA','KA_IB','KAHP','KAHP_DEEPPYR','KAHP_SLOWER','KC','KC_FAST','KDR','KDR_FS','KM','NaF','NaF2','NaF_TCR','NaP','NaPF','NaPF_SS','NaPF_TCR']

# 
# config.py ends here
