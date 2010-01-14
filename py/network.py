# network.py --- 
# 
# Filename: network.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
# Last-Updated: Thu Jan 14 00:31:05 2010 (+0530)
#           By: subhasis ray
#     Update #: 49
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This is for setting up the network. The connectivity is
# stored as a map of maps.
# 
# 
# 
# 

# Change log:
# 
# 
# 
# 

# Code:

from collections import defaultdict
import csv 

class NestedDict(defaultdict):
    """This is a very inefficient but quick-to-code representation for the connectivity matrix."""
    def __init__(self):
	self.default_factory = type(self)


def load_connmatrix(filename):
    """Returns the connectivity matrix as a nested dictionary. The
    file specified by filname is assumed to have first row as the
    column names.

    The entries are like:
    {presynaptic_cell_type_1: {
	postsynaptic_type_1: # of presynaptic_cell_type_1 per postsynaptic_type_1,
	postsynaptic_type_2: # of presynaptic_cell_type_1 per postsynaptic_type_2,
	...},
     presynaptic_cell_type_2:{
	postsynaptic_type_1: # of presynaptic_cell_type_2 per postsynaptic_type_1,
	postsynaptic_type_2: # of presynaptic_cell_type_2 per postsynaptic_type_2,
	...},
    }

    Thus,
    connmatrix = load_connmatrix('file.csv')
    num = connmatrix['presynaptic_cell_type2']['presynaptic_cell_type_2]

    will set the value of num to the number of presynaptic_cell_type_2
    projecting to each postsynaptic_type_2.
    """
    connmatrix = NestedDict()
    reader = csv.reader(file(filename))
    header = reader.next()
    row = 0
    for line in reader:
	if len(line) <= 0:
	    continue
	pre = header[row]
	row += 1
	col = 0
	for entry in line:
	    post = header[col]
	    value = int(entry)
	    connmatrix[pre][post] = value
	    col += 1

    return connmatrix


# 
# network.py ends here
