#!/usr/bin/env python
# plot.py --- 
# 
# Filename: plot.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 24 15:51:57 2009 (+0530)
# Version: 
# Last-Updated: Mon Jun  1 20:44:09 2009 (+0530)
#           By: subhasis ray
#     Update #: 236
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
#  The scaling of MOOSE time must be done with reference to simdt in config.py
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
import os
import pylab
import datetime

working_dir = '/home/subha/src/sim/cortical/'
os.chdir(working_dir)
mus_data_dir = 'py/data/'+datetime.datetime.now().strftime('%Y_%m_%d') + '/'
nrn_data_dir = 'nrn/mydata/'

nrn_data = pylab.loadtxt(nrn_data_dir + 'Vm_ss.plot')
nrn_Vm = nrn_data[:, 1]
nrn_t = nrn_data[:, 0]
mus_Vm = pylab.loadtxt(mus_data_dir + "Vm_ss.plot") * 1e3
mus_t = pylab.linspace(0, len(mus_Vm) * 1e-2, len(mus_Vm))

# for ii in range(1,60):    
#     mus_Vm = pylab.loadtxt(mus_data_dir + "/Vm_comp_" + str(ii) + ".plot") * 1e3
#     mus_t = pylab.linspace(0, len(mus_Vm) * 1e-3, len(mus_Vm))
#     nrn_data = pylab.loadtxt(nrn_data_dir + "/Vm_comp_" + str(ii) + ".plot")
#     nrn_t = nrn_data[:, 0]
#     nrn_Vm = nrn_data[:, 1]
#     pylab.title("Vm_comp_" + str(ii))
#     pylab.plot(nrn_t, nrn_Vm, 'r_', label='nrn_Vm')
#     pylab.plot(mus_t, mus_Vm, 'g.-', label='mus_Vm')
#     pylab.legend()
#     pylab.show()

# from scipy.interpolate import splrep, splev
# smoothness = 3
# order = 2
# # find the knot points
# tckp = splrep( mus_t, mus_Ca, s=3,k=3)

# # evaluate spline, including interpolated points
# mus_Ca_new = splev(nrn_t,tckp)
# print len(mus_Ca_new), len(nrn_Ca)
# # pylab.plot(mus_Ca_new - nrn_Ca)
# # pylab.show()

# # pylab.legend()
# # pylab.show()

# pylab.subplot(2, 1, 1)
pylab.plot(nrn_t, nrn_Vm, 'r-', label='nrn_Vm')
pylab.plot(mus_t, mus_Vm, 'g--', label='mus_Vm')
pylab.legend()
# pylab.subplot(2,1,2)
# pylab.plot(nrn_t, nrn_Ca, label='nrn')
# pylab.plot(mus_t, mus_Ca, label='mus')
# pylab.legend()
# pylab.subplot(2,1,2)
# ratio = mus_gk[2:]*1e-4 / nrn_gk[2:]
# pylab.plot(nrn_t[2:], ratio, 'bx', label='Gk_mus / Gk_nrn')
# pylab.plot(nrn_t, nrn_gk, label='nrn')
# pylab.plot(mus_t, mus_gk, label='mus')
# pylab.legend()
pylab.show()


# 
# plot.py ends here
