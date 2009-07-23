# compare.py --- 
# 
# Filename: compare.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Jul 17 18:01:06 2009 (+0530)
# Version: 
# Last-Updated: Thu Jul 23 23:38:46 2009 (+0530)
#           By: subhasis ray
#     Update #: 116
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

def comparecsv(left, right, epsilon=1e3, row_header=True, col_header=True):
    ret = True
    left_file = open(left, 'rb')
    right_file = open(right, 'rb')
    left_reader = csv.reader(left_file, delimiter=',')
    right_reader = csv.reader(right_file, delimiter=',')
    index = 0
    if col_header:
        left_header = left_reader.next()
        right_header = right_reader.next()
        index = 1

    while True:
        left_row = left_reader.next()
        right_row = right_reader.next()
        if left_row is None and not (right_row is None):
            print left, 'run out of line after', index, 'rows'
            ret = False
            break
        if right_row is None and not (left_row is None):
            print right, 'run out of line after', index, 'rows'
            ret = False
            break
        if len(left_row) != len(right_row):
            print 'No. of columns differ: left - ', len(left_row), 'right -', len(right_row)
            ret = False
            break
        start_col = 0
        if row_header:
            start_col = 1
        col_no = start_col
        for left_col, right_col, in zip(left_row[start_col:], right_row[start_col:]):
            if not almost_equal(float(left_col), float(right_col), epsilon):
                if not row_header:
                    row = str(index)
                else:
                    row = str(left_row[0])
                if not col_header:
                    col = str(col_no)
                else:
                    col = str(left_header[col_no])
                print 'Mismatch in row:%s, column:%s. Values: %f <> %f' % (row, col, left_col, right_col)
                ret = False
        return ret


def compare_tcr():
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

if __name__ == "__main__":
    print comparecsv('py/tcr.txt', 'nrn/TCR')
# 
# compare.py ends here
