# test_hdf5.py --- 
# 
# Filename: test_hdf5.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Oct 23 01:18:37 2010 (+0530)
# Version: 
# Last-Updated: Sat Oct 23 01:39:21 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 42
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

import numpy
import tables

def write_bz2(data):
    bz2_filter = tables.Filters(complib='bzip2', complevel=9)
    bz2_out = tables.openFile('bzip2.h5', mode='w')
    atom = tables.Float64Atom()
    ca = bz2_out.createCArray(bz2_out.root, 'ca', atom, data.shape, filters=bz2_filter)
    ca[:] = data
    bz2_out.close()

def write_zlib(data):
    zlib_filter = tables.Filters(complib='zlib', complevel=9)
    zlib_out = tables.openFile('zlib.h5', mode='w')
    atom = tables.Float64Atom()
    ca = zlib_out.createCArray(zlib_out.root, 'ca', atom, data.shape, filters=zlib_filter)
    ca[:] = data
    zlib_out.close()

def write_lzo(data):
    lzo_filter = tables.Filters(complib='lzo', complevel=9)
    lzo_out = tables.openFile('lzo.h5', mode='w')
    atom = tables.Float64Atom()
    ca = lzo_out.createCArray(lzo_out.root, 'ca', atom, data.shape, filters=lzo_filter)
    ca[:] = data
    lzo_out.close()

def write_uncompressed(data):
    uncompressed_out = tables.openFile('uncompressed.h5', mode='w')
    atom = tables.Float64Atom()
    ca = uncompressed_out.createArray(uncompressed_out.root, 'ca', data)
    ca[:] = data
    uncompressed_out.close()
    
    
def test_compression(size):
    data = numpy.zeros(size)
    data[:size/4] = numpy.random.random_sample(size/4)
    data[size/4:size/2] = numpy.linspace(0, 1, size/4)
    write_bz2(data)
    write_zlib(data)
    write_lzo(data)
    write_uncompressed(data)
    

if __name__ == '__main__':
    test_compression(1000000)


# 
# test_hdf5.py ends here
