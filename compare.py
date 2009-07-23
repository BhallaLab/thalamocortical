# compare.py --- 
# 
# Filename: compare.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Jul 17 18:01:06 2009 (+0530)
# Version: 
# Last-Updated: Tue Jul 21 15:03:35 2009 (+0530)
#           By: subhasis ray
#     Update #: 69
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: compare the properties in neuron and pymoose models printed in a csv file
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

import csv
from numpy import *
def almost_equal(left, right, epsilon=1e-6):
    """check if two floats are almost equal"""
    if left == right:
        return True
    if abs(left) > abs(right):
        return (1 - right / left) < epsilon
    else:
        return ( 1 - left / right) < epsilon
#!almost_equal

if __name__ == '__main__':
    moose_file = open('py/tcr.txt', 'r')
    line = moose_file.readline()
    moose_header = line.split(',')
    nrn_file = open('nrn/tcr.txt', 'r')
    line = nrn_file.readline()
    nrn_header = line.split(',')
    moose_data = loadtxt('py/tcr.txt', delimiter=',', skiprows=1)
    nrn_data = loadtxt('nrn/tcr.txt', delimiter=',', skiprows=1)
    all_matching = True
    for row in range(len(moose_data)):
        for column in range(len(moose_data[row])):
            if not almost_equal(moose_data[row][column], nrn_data[row][column], epsilon=1e-3):
                print 'not matching: comp#', int(moose_data[row][0]), '(', moose_header[column], ')', moose_data[row][column], nrn_data[row][column]
                all_matching = False
    if all_matching: print 'All matching ...'

# 
# compare.py ends here
