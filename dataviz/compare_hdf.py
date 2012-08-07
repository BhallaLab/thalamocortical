# compare_hdf.py --- 
# 
# Filename: compare_hdf.py
# Description: 
# Author: 
# Maintainer: 
# Created: Sun Aug  5 14:07:36 2012 (+0530)
# Version: 
# Last-Updated: Sun Aug  5 16:51:53 2012 (+0530)
#           By: subha
#     Update #: 108
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# compare_hdf.py file1 file2 node
# 
# Recursively compare two hdf5 files under specified node.
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

import h5py as h5
import sys

def compare_recursive(left, right, indent=''):
    ret = 1
    tl = type(left)
    tr = type(right)
    if tl != tr:
        print indent, left, '(', tl,') <>', right, '(', tr, ')'
        return 0
    if isinstance(left, h5.Group):
        lkeys = set(left.keys())
        rkeys = set(right.keys())        
        ldiff = lkeys - rkeys
        if len(ldiff) > 0:
            ret = 0
        for key in ldiff:
            print indent, 'Only on left:', key
        rdiff = rkeys - lkeys
        if len(ldiff) > 0:
            ret = 0
        for key in rdiff:
            print indent, 'Only on right:', key
        inter = lkeys & rkeys
        for key in inter:
            print indent, 'Comparing', key
            ret = compare_recursive(left[key], right[key], indent+' ') and ret
    elif left.attrs['CLASS'] == 'TABLE':
        lval = dict(left.value)
        rval = dict(right.value)
        lkeys = set(lval.keys())
        rkeys = set(rval.keys())
        ldiff = lkeys - rkeys
        if len(ldiff) > 0:
            ret = 0
        for key in ldiff:
            print indent, 'Only on left:', key
        rdiff = rkeys - lkeys
        if len(ldiff) > 0:
            ret = 0
        for key in rdiff:
            print indent, 'Only on right:', key
        inter = lkeys & rkeys
        for key in inter:
            if lval[key] != rval[key]:
                print indent, key,':', lval[key], '<>', rval[key]
                ret = 0
    elif left.attrs['CLASS'] == 'CARRAY':
        lval = left.value
        rval = right.value
        for ii in range(min(len(lval), len(rval))):
            if lval[ii] != rval[ii]:
                print indent, 'row:',ii, ':', lval[ii], '<>', rval[ii]
                ret = 0
        if len(lval) < len(rval):
            ret = 0
            print indent, 'Only on right' 
            for ii in range(len(lval), len(rval)):
                print rval[ii]
        if len(rval) < len(lval):            
            ret = 0
            print indent, 'Only on left'
            for ii in range(len(rval), len(lval)):
                print lval[ii]
    else:
        print indent, 'Could not compare:', left.attrs['CLASS']
    return ret
        

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: %s file1 file2 node\n' % (sys.argv[0])
        print 'compare file1 with file2 from node recursively.'
        sys.exit(0)
    f1 = h5.File(sys.argv[1], 'r')
    f2 = h5.File(sys.argv[2], 'r')
    node = sys.argv[3]
    n1 = f1[node]
    n2 = f2[node]
    print compare_recursive(n1, n2)
    f1.close()
    f2.close()


# 
# compare_hdf.py ends here
