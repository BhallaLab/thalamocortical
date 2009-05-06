#!/usr/bin/env python
# plot.py --- 
# 
# Filename: plot.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 24 15:51:57 2009 (+0530)
# Version: 
# Last-Updated: Wed May  6 18:19:23 2009 (+0530)
#           By: subhasis ray
#     Update #: 182
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
import os
import pylab
import datetime

working_dir = '/home/subha/src/sim/cortical/'
os.chdir(working_dir)
mus_data_dir = 'py/data/'+datetime.datetime.now().strftime('%Y_%m_%d') + '/'
#mus_data_dir = 'py/data/2009_04_24/'
mus_Vm = pylab.loadtxt(mus_data_dir + 'Vm.plot') * 1e3
# mus_Ca = pylab.loadtxt(mus_data_dir + 'Ca.plot') * 1e3
# mus_m = pylab.loadtxt(mus_data_dir + 'm_kahp.plot')
mus_t = pylab.linspace(0, 50, len(mus_Vm))

indices = pylab.array(range(len(mus_Vm)), dtype=int)
mus_indices = indices
nrn_indices = indices * 10
nrn_data_dir = 'nrn/mydata/'

nrn_data = pylab.loadtxt(nrn_data_dir + 'Vm.plot')
# nrn_Ca = pylab.loadtxt(nrn_data_dir + 'Ca.plot')[:, 1]
nrn_Vm = nrn_data[:, 1]
# nrn_m = nrn_data[:, 2]
nrn_t = nrn_data[:, 0]

mus_gk = pylab.loadtxt(mus_data_dir + 'Gk_NaF2.plot')
nrn_gk = nrn_data[:, 2]

# pylab.savetxt('nrn_m.txt', pylab.transpose(( indices, nrn_m[nrn_indices])))
# pylab.savetxt('mus_m.txt', pylab.transpose((mus_indices, mus_m/9.42e-6)))
# pylab.plot(nrn_m[nrn_indices[3:]] , (mus_m[mus_indices[3:]]/9.42e-6))
# pylab.plot( nrn_t, 1/(nrn_m/nrn_Ca), 'bx', label='nrn')
# pylab.plot(mus_t, 1/((mus_m/9.42e-10)/mus_Ca), 'r+', label='mus')

# pylab.plot(mus_t, mus_Ca, 'gx')
# pylab.plot(nrn_t, nrn_Ca*1e2, 'r-')

# pylab.plot(mus_Ca, mus_m/9.42e-6, label='mus')
# pylab.plot(nrn_Ca, nrn_m, label='nrn')
# pylab.plot(nrn_Ca, nrn_t, label='t')
# pylab.plot(mus_t, mus_m/9.42e-6, label='mus')
# pylab.plot(nrn_t, nrn_m, label='nrn')
#pylab.plot(mus_indices[3:], pylab.log(nrn_m[nrn_indices[3:]] / mus_m[mus_indices[3:]]))
#-----
# pylab.plot(mus_t[mus_indices[10:]], pylab.log(mus_m[mus_indices[10:]]), label='mus')
# pylab.plot(nrn_t[nrn_indices[10:]], pylab.log(nrn_m[nrn_indices[10:]]), label='nrn')
#----

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
pylab.plot(nrn_t, nrn_Vm * 1e-2, 'r-', label='nrn_Vm')
# pylab.plot(mus_t, mus_Vm, 'r+', label='mus')
# pylab.legend()
# pylab.subplot(2,1,2)
# pylab.plot(nrn_t, nrn_Ca, label='nrn')
# pylab.plot(mus_t, mus_Ca, label='mus')
# pylab.legend()
# pylab.subplot(2,1,2)
ratio = mus_gk[2:]*1e-4 / nrn_gk[2:]
pylab.plot(nrn_t[2:], ratio, 'bx', label='Gk_mus / Gk_nrn')
# pylab.plot(nrn_t, nrn_gk, label='nrn')
# pylab.plot(mus_t, mus_gk, label='mus')
# pylab.legend()
pylab.show()


# 
# plot.py ends here
