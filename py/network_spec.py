# network_spec.py --- 
# 
# Filename: network_spec.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
# Last-Updated: Mon Mar 29 12:13:40 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 341
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
# 2010-03-29 12:12:22 (+0530) -- moved network.py to network_spec.py
#                                this is currently a backup file.

# Code:

from collections import defaultdict
import csv 
import numpy
from numpy import random

from allowedcomp import ALLOWED_COMP

# may not need nested dict - a 2D dict should be fine: defaultdict(dict)
class NestedDict(defaultdict):
    """This is a very inefficient but quick-to-code representation for
    the connectivity matrix. Good enough for a 14x14 matrix."""
    def __init__(self):
	self.default_factory = type(self)

CELLCOUNT = {
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


def load_connmap(filename):
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
    connmap = load_connmap('file.csv')
    num = connmap['presynaptic_cell_type2']['presynaptic_cell_type_2]

    will set the value of num to the number of presynaptic_cell_type_2
    projecting to each postsynaptic_type_2.
    """
    connmap = defaultdict(dict)
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
	    connmap[pre][post] = value
	    col += 1

    return connmap

def test_load_connmap(filename='connmatrix.txt'):
    matrix = load_connmap(filename)
    for key, value in matrix.items():
        print '\n#', key, '#'
        for key1, value1 in value.items():
            print '\t', key1, '\t:', value1


CONN_MATRIX = load_connmap('connmatrix.txt')

NUM_ALLOWED_COMP = defaultdict(dict)
NUM_ALLOWED_COMP['SupPyrRS']['SupPyrRS'] = 36
NUM_ALLOWED_COMP['SupPyrRS']['SupPyrFRB'] = 36
NUM_ALLOWED_COMP['SupPyrRS']['SupBasket']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['SupAxoaxonic']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['SupLTS']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['SpinyStellate']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['TuftedIB']  =  8
NUM_ALLOWED_COMP['SupPyrRS']['TuftedRS']  =  8
NUM_ALLOWED_COMP['SupPyrRS']['DeepBasket']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['DeepAxoaxonic']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['DeepLTS']  = 24
NUM_ALLOWED_COMP['SupPyrRS']['NontuftedRS']  =  7

NUM_ALLOWED_COMP['SupPyrFRB']['SupPyrRS'] = 36
NUM_ALLOWED_COMP['SupPyrFRB']['SupPyrFRB'] = 36
NUM_ALLOWED_COMP['SupPyrFRB']['SupBasket'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['SupAxoaxonic'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['SupLTS'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['SpinyStellate'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['TuftedIB'] =  8
NUM_ALLOWED_COMP['SupPyrFRB']['TuftedRS'] =  8
NUM_ALLOWED_COMP['SupPyrFRB']['DeepBasket'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['DeepLTS'] = 24
NUM_ALLOWED_COMP['SupPyrFRB']['NontuftedRS'] =  7

NUM_ALLOWED_COMP['SupBasket']['SupPyrRS']   = 11
NUM_ALLOWED_COMP['SupBasket']['SupPyrFRB']   = 11
NUM_ALLOWED_COMP['SupBasket']['SupBasket']     = 24
NUM_ALLOWED_COMP['SupBasket']['SupAxoaxonic']     = 24
NUM_ALLOWED_COMP['SupBasket']['SupLTS']      = 24
NUM_ALLOWED_COMP['SupBasket']['SpinyStellate']   =  5

NUM_ALLOWED_COMP['SupLTS']['SupPyrRS'] = 53
NUM_ALLOWED_COMP['SupLTS']['SupPyrFRB'] = 53
NUM_ALLOWED_COMP['SupLTS']['SupBasket'] = 40
NUM_ALLOWED_COMP['SupLTS']['SupAxoaxonic'] = 40
NUM_ALLOWED_COMP['SupLTS']['SupLTS'] = 40
NUM_ALLOWED_COMP['SupLTS']['SpinyStellate'] = 40
NUM_ALLOWED_COMP['SupLTS']['TuftedIB'] = 40
NUM_ALLOWED_COMP['SupLTS']['TuftedRS'] = 40
NUM_ALLOWED_COMP['SupLTS']['DeepBasket'] = 20
NUM_ALLOWED_COMP['SupLTS']['DeepAxoaxonic'] = 20
NUM_ALLOWED_COMP['SupLTS']['DeepLTS'] = 20
NUM_ALLOWED_COMP['SupLTS']['NontuftedRS'] = 29

NUM_ALLOWED_COMP['SpinyStellate']['SupPyrRS'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['SupPyrFRB'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['SupBasket'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['SupAxoaxonic'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['SupLTS'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['SpinyStellate'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['TuftedIB'] = 12
NUM_ALLOWED_COMP['SpinyStellate']['TuftedRS'] = 12
NUM_ALLOWED_COMP['SpinyStellate']['DeepBasket'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['DeepLTS'] = 24
NUM_ALLOWED_COMP['SpinyStellate']['NontuftedRS'] =  5

NUM_ALLOWED_COMP['TuftedIB']['SupPyrRS'] = 13
NUM_ALLOWED_COMP['TuftedIB']['SupPyrFRB'] = 13
NUM_ALLOWED_COMP['TuftedIB']['SupBasket'] = 24
NUM_ALLOWED_COMP['TuftedIB']['SupAxoaxonic'] = 24
NUM_ALLOWED_COMP['TuftedIB']['SupLTS'] = 24
NUM_ALLOWED_COMP['TuftedIB']['SpinyStellate'] = 24
NUM_ALLOWED_COMP['TuftedIB']['TuftedIB'] = 46
NUM_ALLOWED_COMP['TuftedIB']['TuftedRS'] = 46
NUM_ALLOWED_COMP['TuftedIB']['DeepBasket'] = 24
NUM_ALLOWED_COMP['TuftedIB']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['TuftedIB']['DeepLTS'] = 24
NUM_ALLOWED_COMP['TuftedIB']['NontuftedRS'] = 43

NUM_ALLOWED_COMP['TuftedRS']['SupPyrRS'] = 13
NUM_ALLOWED_COMP['TuftedRS']['SupPyrFRB'] = 13
NUM_ALLOWED_COMP['TuftedRS']['SupBasket'] = 24
NUM_ALLOWED_COMP['TuftedRS']['SupAxoaxonic'] = 24
NUM_ALLOWED_COMP['TuftedRS']['SupLTS'] = 24
NUM_ALLOWED_COMP['TuftedRS']['SpinyStellate'] = 24
NUM_ALLOWED_COMP['TuftedRS']['TuftedIB'] = 46
NUM_ALLOWED_COMP['TuftedRS']['TuftedRS'] = 46
NUM_ALLOWED_COMP['TuftedRS']['DeepBasket'] = 24
NUM_ALLOWED_COMP['TuftedRS']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['TuftedRS']['DeepLTS'] = 24
NUM_ALLOWED_COMP['TuftedRS']['NontuftedRS'] = 43

NUM_ALLOWED_COMP['DeepBasket']['SpinyStellate'] = 5 
NUM_ALLOWED_COMP['DeepBasket']['TuftedIB'] =  8
NUM_ALLOWED_COMP['DeepBasket']['TuftedRS'] =  8
NUM_ALLOWED_COMP['DeepBasket']['DeepBasket'] = 24
NUM_ALLOWED_COMP['DeepBasket']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['DeepBasket']['DeepLTS'] = 24
NUM_ALLOWED_COMP['DeepBasket']['NontuftedRS'] =  8

NUM_ALLOWED_COMP['DeepLTS']['SupPyrRS'] = 53
NUM_ALLOWED_COMP['DeepLTS']['SupPyrFRB'] = 53
NUM_ALLOWED_COMP['DeepLTS']['SupBasket'] = 20
NUM_ALLOWED_COMP['DeepLTS']['SupAxoaxonic'] = 20
NUM_ALLOWED_COMP['DeepLTS']['SupLTS'] = 20
NUM_ALLOWED_COMP['DeepLTS']['SpinyStellate'] = 40
NUM_ALLOWED_COMP['DeepLTS']['TuftedIB'] = 40
NUM_ALLOWED_COMP['DeepLTS']['TuftedRS'] = 40
NUM_ALLOWED_COMP['DeepLTS']['DeepBasket'] = 40
NUM_ALLOWED_COMP['DeepLTS']['DeepAxoaxonic'] = 40
NUM_ALLOWED_COMP['DeepLTS']['DeepLTS'] = 40
NUM_ALLOWED_COMP['DeepLTS']['NontuftedRS'] = 29


NUM_ALLOWED_COMP['TCR']['SupPyrRS'] = 24
NUM_ALLOWED_COMP['TCR']['SupPyrFRB'] = 24
NUM_ALLOWED_COMP['TCR']['SupBasket'] = 12
NUM_ALLOWED_COMP['TCR']['SupAxoaxonic'] = 12
NUM_ALLOWED_COMP['TCR']['SpinyStellate'] = 52
NUM_ALLOWED_COMP['TCR']['TuftedIB'] =  9
NUM_ALLOWED_COMP['TCR']['TuftedRS'] =  9
NUM_ALLOWED_COMP['TCR']['DeepBasket'] = 12
NUM_ALLOWED_COMP['TCR']['DeepAxoaxonic'] = 12
NUM_ALLOWED_COMP['TCR']['nRT'] = 12
NUM_ALLOWED_COMP['TCR']['NontuftedRS'] =  5

NUM_ALLOWED_COMP['nRT']['TCR'] = 11
NUM_ALLOWED_COMP['nRT']['nRT'] = 53

NUM_ALLOWED_COMP['NontuftedRS']['SupPyrRS'] = 4
NUM_ALLOWED_COMP['NontuftedRS']['SupPyrFRB'] =  4
NUM_ALLOWED_COMP['NontuftedRS']['SupBasket'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['SupAxoaxonic'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['SupLTS'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['SpinyStellate'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['TuftedIB'] = 46
NUM_ALLOWED_COMP['NontuftedRS']['TuftedRS'] = 46
NUM_ALLOWED_COMP['NontuftedRS']['DeepBasket'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['DeepAxoaxonic'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['DeepLTS'] = 24
NUM_ALLOWED_COMP['NontuftedRS']['TCR'] = 90
NUM_ALLOWED_COMP['NontuftedRS']['nRT'] = 12
NUM_ALLOWED_COMP['NontuftedRS']['NontuftedRS'] = 43

import moose
class Population(moose.Neutral):
    """Cell populations - handles setting up connections between
    populations based on connection matrix."""
    def __init__(self, path, cellclass, number, prefix):
        """Generate 'number' cells of class 'cellclass' under the
        neutral object at 'path', the n-th cell is called
        {prefix}_n"""
        moose.Neutral.__init__(self, path)
        self.cellList = []
        self.cellType = cellclass.__name__
        if not isinstance(cellclass, type):
            raise Error, 'Population.__init__(): Please provide a class object for crllclass'        
        for count in range(number):
            cell_name = prefix + '_' + str(count)
            cell_instance = cellclass(cell_name, self)
            self.cell.append(cell_instance)


    def populationConnect(self, post):
        """Create a connection between this population and the
        population specified by 'post'"""
        num_pre_per_post = CONNMAP[self.cellType][post.cellType]
        # This gets a 2_D matrix whose row[i] is the array of indices
        # of the pre-synaptic cells for i-th post-synaptic cell.
        precell_indices = numpy.random.randint(0, high=len(self.cellList), size=(len(post.cellList), num_pre_per_post))

        allowed_comp_list = array(ALLOWED_COMP[self.cellType][post.cellType])
        target_comp_indices = numpy.random.randint(0, high=len(allowed_comp_list), size=(len(post.cellList), len(num_pre_per_post)))
        for ii in range(len(post.cellList)):
            target_comp_list = allowed_comp_list[target_comp_indices[ii]] # list containing the target compartment for each presynaptic cell
            for jj in range(len(num_pre_per_post)):
                precell_index = precell_indices[ii][jj]
                precell = self.cellList[precell_index]
                postcell = post.cellList[ii]
                post_syn_comp_index = target_comp_indices[jj]
                post_syn_comp = postcell.comp[post_syn_comp_index]
                pre_syn_comp = precell.comp[precell.presyn]
            


def get_synaptic_compmap(pre_cell_per_post_cell, allowed_comps):
    '''returns a list containing the post-synaptic compartment no. for
    each presynaptic cell.

    pre_cell_per_post_cell is the number of presynaptic cells
    connecting to each post-synaptic cell.

    allowed_comps is a list containing the indices of the compartments
    on the post synaptic cell that are valid targets of a synapse.
    '''
    comp_map = numpy.random.randint(0, high=len(allowed_comps), size=len(pre_cell_per_post_cell))
    # for jj in range(len(pre_cell_per_post_cell)):
    #     target_allowed_comp_index = random.randint(0, len(allowed_comps) - 1)
    #     comp_map.append(allowed_comps[jj])
    return comp_map


def test_num_allowedcom():
    for pre, dest in ALLOWED_COMP.items():
        for post, comps in dest.items():
            if len(comps) != NUM_ALLOWED_COMP[pre][post]:
                print 'Error:', pre, post


if __name__ == '__main__':
    # test_load_connmap()
    test_num_allowedcom()

                             
                             
# 
# network_spec.py ends here
