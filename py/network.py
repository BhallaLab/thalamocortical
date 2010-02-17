# network.py --- 
# 
# Filename: network.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
# Last-Updated: Wed Feb 17 16:53:34 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 165
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
# 2010-01-13 initial version to load csv
#
# 2010-01-14 added cell count for the original model
# 
# 

# Code:

from collections import defaultdict
import csv 
import numpy

# may not need nested dict - a 2D dict should be fine: defaultdict(dict)
class NestedDict(defaultdict):
    """This is a very inefficient but quick-to-code representation for
    the connectivity matrix. Good enough for a 14x14 matrix."""
    def __init__(self):
	self.default_factory = type(self)

cellcount = {
    'SupPyrRS': 1000,
    'SupPyrFRB': 50,
    'SupBasket': 90,
    'SupAxoaxonic': 90,
    'SupLTS': 90,
    'SpinyStellate': 240,
    'TuftedIB': 800,
    'TuftedRS': 200,
    'NontuftedRS': 500,
    'DeepBasket': 100,
    'DeepAxoaxonic': 100,
    'DeepLTS': 100,
    'TCR': 100,
    'nRT': 100
}

class ConnectionData:
    def __init__(self):
        self.celltype = ["SupPyrRS",
                         "SupPyrFRB",
                         "SupBasket",
                         "SupAxoaxonic",
                         "SupLTS",
                         "SpinyStellate",
                         "TuftedIB",
                         "TuftedRS",
                         "DeepBasket",
                         "DeepAxoaxonic",
                         "DeepLTS",
                         "NontuftedRS",
                         "TCR",
                         "nRT"]
        
        # Conn matrix is the representation of the square matrix whose
        # element [i][j] represents how many celltype[i] connect to a
        # single cell of celltype[j]
        self.matrix = numpy.array[[50, 50, 90, 90, 90,  3, 60, 60, 30, 30, 30,  3,  0,  0],
                                  [ 5,  5,  5,  5,  5,  1,  3,  3,  3,  3,  3,  1,  0,  0],
                                  [20, 20, 20, 20, 20, 20,  0,  0,  0,  0,  0,  0,  0,  0],
                                  [20, 20,  0,  0,  0,  5,  5,  5,  0,  0,  0,  5,  0,  0],
                                  [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,  0,  0],
                                  [20, 20, 20, 20, 20, 30, 20, 20, 20, 20, 20, 20,  0,  0],
                                  [ 2,  2, 20, 20, 20, 20, 50, 20, 20, 20, 20, 20,  0,  0],
                                  [ 2,  2, 20, 20, 20, 20, 20, 10, 20, 20, 20, 20,  0,  0],
                                  [ 0,  0,  0,  0,  0, 20, 20, 20, 20, 20, 20, 20,  0,  0],
                                  [ 5,  5,  0,  0,  0,  5,  5,  5,  0,  0,  0,  5,  0,  0],
                                  [10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20,  0,  0],
                                  [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20],
                                  [10, 10, 10, 10,  0,  0, 10, 10, 20, 10,  0, 10,  0, 25],
                                  [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 15, 10]]

    def precell_count(self, precelltype, postcelltype):
        """This returns the number of cells of precelltype connect to
        one cell of postcelltype."""
        pre_index = self.matrix.index(precelltype)
        post_index = self.matrix.index(postcelltype)
        return self.matrix[pre_index][post_index]
# Till here - ConnectionData may be unnecessary if we are storing it in file                    


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
#    connmatrix = NestedDict()
    connmatrix = defaultdict(dict)
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

def test_load_connmatrix(filename='connmatrix.txt'):
    matrix = load_connmatrix(filename)
    for key, value in matrix.items():
        print '\n#', key, '#'
        for key1, value1 in value.items():
            print '\t', key1, '\t:', value1

import moose
class Population(moose.Neutral):
    """Cell populations - handles setting up connections between
    populations based on connection matrix."""
    def __init__(self, path, cellclass, number, prefix):
        """Generate 'number' cells of class 'cellclass' under the
        neutral object at 'path', the n-th cell is called
        {prefix}_n"""
        moose.Neutral.__init__(self, path)
        self.cell = []
        for count in range(number):
            self.cell.append(moose.Cell(prefix + '_' + str(count), self))


    def populationConnect(self, post, pre_count):
        """Create a connection between this population and the
        population specified by 'post' where 'pre_count' cells
        from this population will be connected to one cell in 'post'
        population"""
        # Select the pre_count cells from this presynaptic population.
        # pre = 
        pass
        
if __name__ == '__main__':
    test_load_connmatrix()

                             
                             
# 
# network.py ends here
