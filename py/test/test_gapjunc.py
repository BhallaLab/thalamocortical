# test_gapjunc.py --- 
# 
# Filename: test_gapjunc.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Oct  5 14:05:08 2010 (+0530)
# Version: 
# Last-Updated: Tue Oct  5 14:09:26 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 18
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
ctx = moose.PyMooseBase.getContext()

cell_a = moose.Compartment('a')
cell_b = moose.Compartment('b')
gj = moose.Compartment('gj')
cell_a.diameter = 1e-6
cell_a.len = 20e-6
cell_b.diameter = 1e-6
cell_b.len = 20e-6
cell_a.Rm = 1e5
cell_a.Ra = 1e6
cell_b.Rm = 1e5
cell_b.Ra = 1e6
gj.Ra = 1e6
cell_a.connect('raxial', gj, 'axial')
gj.connect('raxial', cell_b, 'axial')
ctx.reset()
ctx.step(1000)



# 
# test_gapjunc.py ends here
