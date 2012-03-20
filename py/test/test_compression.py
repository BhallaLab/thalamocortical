# test_compression.py --- 
# 
# Filename: test_compression.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Mon Oct 25 14:52:02 2010 (+0530)
# Version: 
# Last-Updated: Mon Oct 25 18:55:40 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 158
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

# Code:

import os
import sys
from datetime import datetime
from collections import defaultdict
import numpy
import tables
import pylab

def save_h5(filepath, complib, diff):
    """Read series of floating point number from filepath and save it
    as hdf5. The extension is .h5 if diff is False. If diff is True
    the file ends with _diff.h5

    filepath -- path of the file to be read.

    complib -- compression library, passed on to PyTables. Can be
    zlib, lzo, gzip, bzip2.

    diff -- save the original numbers or the first differences.
    """
    data = numpy.loadtxt(filepath, dtype=float)
    if diff:
        new_file_path = '%s_diff_%s.h5' % (filepath, complib)
        to_save = numpy.diff(data)        
    else:
        new_file_path = '%s_%s.h5' % (filepath, complib)
        to_save = data
    fltr = tables.Filters(complevel=9, complib=complib)
    new_file = tables.openFile(new_file_path, mode='w')
    carray = new_file.createCArray(new_file.root, 'array', tables.FloatAtom(), numpy.shape(to_save), filters=fltr)
    if diff:
        carray.start = float(data[0])
    carray[:] = to_save
    new_file.close()
    return new_file_path

def test_compression(directory):
    complibs = ['bzip2', 'zlib', 'lzo']
    file_list = defaultdict(list)
    time_list = defaultdict(list)
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        print 'Processing', filepath
        for complib in complibs:
            start = datetime.now()
            new_file_path = save_h5(filepath, complib, False)
            end = datetime.now()
            delta = end - start
            file_list[complib].append(new_file_path)
            time_list[complib].append(delta.seconds + delta.microseconds * 1e-6)

            start = datetime.now()
            new_file_path = save_h5(filepath, complib, True)
            end = datetime.now()
            delta = end - start
            file_list[complib+'_diff'].append(new_file_path)
            time_list[complib+'_diff'].append(delta.seconds + delta.microseconds * 1e-6)

    return (file_list, time_list)
            
def plot_size_and_time(file_list, time_list):
    sizes = []
    times = []
    labels = []
    for key in file_list.keys():
        total_size = 0
        total_time = 0
        for f in file_list[key]:
            total_size += os.path.getsize(f)        
        for t in time_list[key]:
            total_time += t
        sizes.append(total_size)
        times.append(total_time)
        labels.append(key)
        
    pylab.subplot(211)
    width = 0.5
    xlocations = numpy.array(range(len(labels))) + width
    pylab.bar(xlocations, sizes, width=width)
    ticklabels = ['%s (%g MB)' % (labels[ii], sizes[ii]/2**20) for ii in range(len(labels))]
    pylab.xticks(xlocations + width/2.0, ticklabels)
    pylab.title('File sizes')

    pylab.subplot(212)
    pylab.bar(xlocations, times, width=width)
    ticklabels = ['%s (%g s)' % (labels[ii], times[ii]) for ii in range(len(labels))]
    pylab.xticks(xlocations + width/2.0, ticklabels)
    pylab.title('Time to process')

    pylab.show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        directory = '.'
    (files, times) = test_compression(directory)
    plot_size_and_time(files, times)

# 
# test_compression.py ends here
