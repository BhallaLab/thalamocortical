#!/usr/bin/env python
# ig_traubnet.py --- 
# 
# Filename: ig_traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Sep 16 16:19:39 2010 (+0530)
# Version: 
# Last-Updated: Thu Oct  7 18:08:27 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 924
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# igraph version of cell-cell graph generation for traub 2005 model
# 
# Initially the celltype graph was built from manually entered data in
# connmatrix.txt. Which contained an NxN matrix of all the cell
# types. The entry[i][j] had the number of presynaptic cells of
# celltype[i] connecting to each postsynaptic cell of celltype[j].  In
# addition, allowedcomp.py contained a dict ALLOWED_COMP which had a
# map between pairs of celltypes to a list of comprtment
# numbers. These are the compartment numbers in the second cell on
# which the first cell of the pair could form synapse.
#

# Change log:
#
# 2010-09-23 11:17:55 (+0530) -- now the celltype_graph has edges
# which have multiple attributes for all types of
# synapse. taugabafast, taugabaslow,
# 
# 2010-10-05 19:18:49 (+0530) -- changed _generate_celltype_graph and
# _generate_cell_graph to public methods generate_celltype_graph and
# generate_cell_graph. The idea is, I may want to manipulate the
# celltype-graph before generating the cell-graph - especially for
# scaling the synaptic conductances.
# 

# Code:

import sys
import os
from datetime import datetime
import numpy
import igraph as ig
import config
# synapse.py is imported for checking against the older version of the synaptic data.
import synapse
# allowedcomp is imported for verification with older data structure only
import allowedcomp

# The cell classes
from spinystellate import SpinyStellate
from suppyrRS import SupPyrRS
from suppyrFRB import SupPyrFRB
from supbasket import SupBasket
from supaxoaxonic import SupAxoaxonic
from supLTS import SupLTS
from tuftedIB import TuftedIB
from deepbasket import DeepBasket
from deepaxoaxonic import DeepAxoaxonic
from deepLTS import DeepLTS
from tuftedRS import TuftedRS
from nontuftedRS import NontuftedRS
from tcr import TCR
from nRT import nRT


class TraubFullNetData(object):
    """Information about connectivity.

    This is hand coded datastructure for the Traub model network
    structure. In absense of the cell-type network graph or other
    files, this will be used as the master data.
    
    `celltype` -- list of celltypes
    
    `cellcount` -- list of cell-counts for each cell
                   type. cellcount[i] is the number of cells of type celltype[i]
                   present in the network.

    `pre_post_ratio` -- a square matrix represented by a list of lists
                        containing number of presynaptic cell of each
                        type connecting to each post-synaptic cell of
                        each type. pre_post_ratio[i][j] is the number
                        of cells of type celltype[i] forming synapse
                        on each cell of type celltype[j].

    """
    def __init__(self):
        self.celltype = ['SupPyrRS',
                         'SupPyrFRB',
                         'SupBasket',
                         'SupAxoaxonic',
                         'SupLTS',
                         'SpinyStellate',                          
                         'TuftedIB',
                         'TuftedRS',
                         'DeepBasket',
                         'DeepAxoaxonic',
                         'DeepLTS',
                         'NontuftedRS',
                         'TCR',
                         'nRT']

        self.cellcount = [1000,
                          50,
                          90,
                          90,
                          90,
                          240,
                          800,
                          200,
                          100,
                          100,
                          100,
                          500,
                          100,
                          100]

        
        self.pre_post_ratio = [[50,     50,     90,     90,     90,     3,      60,     60,     30,     30,     30,     3,      0,      0],
                               [5,      5,      5,	5,	5,	1,	3,	3,	3,	3,	3,	1,	0,	0],
                               [20,	20,	20,	20,	20,	20,	0,	0,	0,	0,	0,	0,	0,	0],
                               [20,	20,	0,	0,	0,	5,	5,	5,	0,	0,	0,	5,	0,	0],
                               [20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	0,	0],
                               [20,	20,	20,	20,	20,	30,	20,	20,	20,	20,	20,	20,	0,	0],
                               [2,	2,	20,	20,	20,	20,	50,	20,	20,	20,	20,	20,	0,	0],
                               [2,	2,	20,	20,	20,	20,	20,	10,	20,	20,	20,	20,	0,	0],
                               [0,	0,	0,	0,	0,	20,	20,	20,	20,	20,	20,	20,	0,	0],
                               [5,	5,	0,	0,	0,	5,	5,	5,	0,	0,	0,	5,	0,	0],
                               [10,	10,	10,	10,	10,	20,	20,	20,	20,	20,	20,	20,	0,	0],
                               [10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	20,	20,	20],
                               [10,	10,	10,	10,	0,	0,	10,	10,	20,	10,	0,	10,	0,	25],
                               [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	15,	10]]

        self.tau_ampa = [[   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3,    0.0,    0.0],
                         [   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3,    0.0,    0.0],
                         [   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3,    0.0,    0.0],
                         [   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0],
                         [   2.0e-3,    2.0e-3,  0.8e-3,  0.8e-3, 1.0e-3, 2.0e-3, 2.0e-3, 2.0e-3,   0.8e-3,   0.8e-3,  1.0e-3,    2.0e-3, 2.0e-3, 2.0e-3],
                         [   2.0e-3,    2.0e-3,  1.0e-3,  1.0e-3,    0.0, 2.0e-3, 2.0e-3, 2.0e-3,   1.0e-3,   1.0e-3,     0.0,    2.0e-3,    0.0, 2.0e-3],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,    0.0]]

        # ta_nmda[suppyrrs][suppyrrs] = 130.0 ms according to paper, but 130.5 ms in code
        self.tau_nmda = [[ 130.5e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
                         [ 130.0e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [   130e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
                         [   130e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
                         [   130e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0],
                         [   130e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3, 100.0e-3,    130e-3, 130e-3, 100e-3],
                         [   130e-3,    130e-3,  100e-3,  100e-3,    0.0, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,      0.0,    130e-3,    0.0, 150e-3],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,      0.0,       0.0,    0.0,    0.0]]

        self.tau_gaba = [[      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [     6e-3,      6e-3,    3e-3,    3e-3,   3e-3,  6e-3,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [     6e-3,      6e-3,     0.0,     0.0,    0.0,  6e-3,   6e-3,   6e-3,      0.0,      0.0,     0.0,      6e-3,    0.0,  0.0 ],
                         [    20e-3,     20e-3,   20e-3,   20e-3,  20e-3, 20e-3,  20e-3,  20e-3,    20e-3,    20e-3,   20e-3,     20e-3,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,  6e-3,   6e-3,   6e-3,     3e-3,     3e-3,    3e-3,      6e-3,    0.0,  0.0 ],
                         [     6e-3,      6e-3,     0.0,     0.0,    0.0,  6e-3,   6e-3,   6e-3,      0.0,      0.0,     0.0,      6e-3,    0.0,  0.0 ],
                         [    20e-3,     20e-3,   20e-3,   20e-3,  20e-3, 20e-3,  20e-3,  20e-3,    20e-3,    20e-3,   20e-3,     20e-3,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,    0.0,  0.0 ],
                         [      0.0,       0.0,     0.0,     0.0,    0.0,   0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0, 3.3e-3, 9e-3 ]]

        # nRT_tau_gaba_slow - first entry nRT -> TCR, nRT->nRT
        self.nRT_TCR_tau_gaba_slow = 10e-3 
        self.nRT_nRT_tau_gaba_slow = 44.5e-3

        self.g_ampa_baseline = [
            [  0.25e-9,   0.25e-9,    3e-9,    3e-9,   2e-9, 0.1e-9, 0.1e-9, 0.1e-9,     1e-9,     1e-9,    1e-9,    0.5e-9,     0.0,     0.0 ],
            [  0.25e-9,   0.25e-9,    3e-9,    3e-9,   2e-9, 0.1e-9, 0.1e-9, 0.1e-9,     1e-9,     1e-9,    1e-9,    0.5e-9,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [     1e-9,      1e-9,    1e-9,    1e-9,   1e-9,   1e-9,   1e-9,   1e-9,     1e-9,     1e-9,    1e-9,      1e-9,     0.0,     0.0 ],
            [   0.5e-9,    0.5e-9,    1e-9,    1e-9,   1e-9, 0.5e-9,   2e-9,   2e-9,     3e-9,     3e-9,    2e-9,      2e-9,     0.0,     0.0 ],
            [   0.5e-9,    0.5e-9,    1e-9,    1e-9,   1e-9,  0.5-9,   1e-9,   1e-9,     3e-9,     3e-9,    2e-9,      1e-9,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ],
            [   0.5e-9,    0.5e-9,    1e-9,    1e-9,   1e-9, 0.5e-9,   1e-9,   1e-9,     3e-9,     3e-9,    2e-9,      1e-9, 0.75e-9,  0.5e-9 ],
            [   0.5e-9,    0.5e-9,  0.1e-9,  0.1e-9,    0.0,   1e-9, 1.5e-9, 1.5e-9,   1.5e-9,     1e-9,     0.0,      1e-9,     0.0, 0.75e-9 ],
            [      0.0,       0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,      0.0,      0.0,     0.0,       0.0,     0.0,     0.0 ]]

        self.g_nmda_baseline = [
            [ 0.025e-9,  0.025e-9, 0.15e-9, 0.15e-9, 0.15e-9, 0.01e-9, 0.01e-9, 0.01e-9,   0.1e-9,   0.1e-9, 0.15e-9,   0.05e-9,      0.0,     0.0 ],
            [ 0.025e-9,  0.025e-9, 0.15e-9, 0.15e-9, 0.15e-9, 0.01e-9, 0.01e-9, 0.01e-9,   0.1e-9,   0.1e-9, 0.15e-9,   0.05e-9,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [   0.1e-9,    0.1e-9, 0.15e-9, 0.15e-9, 0.15e-9,  0.1e-9,  0.1e-9,  0.1e-9,  0.15e-9,  0.15e-9, 0.15e-9,    0.1e-9,      0.0,     0.0 ],
            [  0.05e-9,   0.05e-9, 0.15e-9, 0.15e-9, 0.15e-9, 0.05e-9,  0.2e-9,  0.2e-9,  0.15e-9,  0.15e-9, 0.15e-9,    0.2e-9,      0.0,     0.0 ],
            [  0.05e-9,   0.05e-9, 0.15e-9, 0.15e-9, 0.15e-9,  0.05-9,  0.1e-9,  0.1e-9,   0.1e-9,   0.1e-9,  0.1e-9,    0.1e-9,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ],
            [  0.05e-9,   0.05e-9,  0.1e-9,  0.1e-9,  0.1e-9, 0.05e-9,  0.1e-9,  0.1e-9,   0.1e-9,   0.1e-9,  0.1e-9,    0.1e-9, 0.075e-9, 0.05e-9 ],
            [  0.05e-9,   0.05e-9, 0.01e-9, 0.01e-9,     0.0,  0.1e-9, 0.15e-9, 0.15e-9,   0.1e-9,   0.1e-9,     0.0,    0.1e-9,      0.0, 0.15e-9 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,      0.0,     0.0 ]]

        # for each nRT->TCR connections, g_gaba_baseline is taken from a uniform distribution in the range 0.7e-9 to 2.1e-9 Siemens.
        self.g_gaba_baseline = [
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [   1.2e-9,    1.2e-9,  0.2e-9,  0.2e-9,  0.5e-9,  0.1e-9,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [   1.2e-9,    1.2e-9,     0.0,     0.0,     0.0,  0.1e-9,    1e-9,    1e-9,      0.0,      0.0,     0.0,      1e-9, 0.0,    0.0 ],
            [  0.01e-9,   0.01e-9, 0.01e-9, 0.01e-9, 0.05e-9, 0.01e-9, 0.02e-9, 0.02e-9,  0.01e-9,  0.01e-9, 0.05e-9,   0.01e-9, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,  1.5e-9,  0.7e-9,  0.7e-9,   0.2e-9,   0.2e-9,  0.7e-9,    0.7e-9, 0.0,    0.0 ],
            [     1e-9,      1e-9,     0.0,     0.0,     0.0,  1.5e-9,    1e-9,    1e-9,      0.0,      0.0,     0.0,      1e-9, 0.0,    0.0 ],
            [  0.01e-9,   0.01e-9, 0.01e-9, 0.01e-9, 0.05e-9, 0.01e-9, 0.05e-9, 0.02e-9,  0.01e-9,  0.01e-9, 0.05e-9,   0.01e-9, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0,    0.0 ],
            [      0.0,       0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0, 0.0, 0.3e-9 ]]

        self.allowed_comps = [
            [# SupPyrRS
                [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26, 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25, 34,35,36,37], # SupPyrRS
                [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26, 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25, 34,35,36,37], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SpinyStellate
                [39,40,41,42,43,44,45,46], # TuftedIB
                [39,40,41,42,43,44,45,46], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # DeepLTS
                [38,39,40,41,42,43,44], # NontuftedRS
                [], # TCR 
                []], # nRT
            [ # SupPyrFRB
                [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26, 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25, 34,35,36,37], # SupPyrRS
                [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26, 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25, 34,35,36,37], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36, 44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SpinyStellate
                [39,40,41,42,43,44,45,46], # TuftedIB
                [39,40,41,42,43,44,45,46], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [38,39,40,41,42,43,44], # NontuftedRS
                [],
                []],
            [# SupBasket
                [1,2,3,4,5,6,7,8,9,38,39], # SupPyrRS
                [1,2,3,4,5,6,7,8,9,38,39], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupLTS
                [1,2,15,28,41], # SpinyStellate
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                []],
            [ # SupAxoaxonic
                [69],
                [69],
                [],
                [],
                [],
                [54],
                [56],
                [56],
                [],
                [],
                [],
                [45],
                [],
                []],
            [ # SupLTS
                [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68] , # SupPyrRS
                [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68], # SupPyrFRB
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # SupBasket
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # SupAxoaxonic
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # SupLTS
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # SpinyStellate
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55], # TuftedIB
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55], # TuftedRS
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # DeepBasket
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # DeepAxoaxonic
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # DeepLTS
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44], # NontuftedRS
                [],
                []],
            [ #SpinyStellate
                [ 2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,26,27,28,29,30,31,32,33], # SupPyrRS
                [ 2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,26,27,28,29,30,31,32,33], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SpinyStellate
                [7,8,9,10,11,12,36,37,38,39,40,41], # TuftedIB
                [7,8,9,10,11,12,36,37,38,39,40,41], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [37,38,39,40,41], # NontuftedRS
                [],
                []],
            [ #TuftedIB
                [40,41,42,43,44,45,46,47,48,49,50,51,52], #SupPyrRS
                [40,41,42,43,44,45,46,47,48,49,50,51,52], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SpinyStellate
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedIB
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44], # NontuftedRS
                [],
                []],
            [ # TuftedRS
                [40,41,42,43,44,45,46,47,48,49,50,51,52], # SupPyrRS
                [40,41,42,43,44,45,46,47,48,49,50,51,52], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SpinyStellate
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedIB
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44], # NontuftedRS
                [],
                []],
            [ # DeepBasket
                [],
                [],
                [],
                [],
                [],
                [1,2,15,28,41], # SpinyStellate
                [1,2,3,4,5,6,35,36], # TuftedIB
                [1,2,3,4,5,6,35,36], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [1,2,3,4,5,6,35,36], # NontuftedRS
                [],
                []],
            [ # DeepAxoaxonic
                [69],
                [69],
                [],
                [],
                [],
                [54],
                [56],
                [56],
                [],
                [],
                [],
                [45],
                [],
                []],
            [ # DeepLTS 
                [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68], # SupPyrRS
                [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68], # SupPyrFRB
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # SupBasket
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # SupAxoaxonic
                [8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,47,48,49,50,51], # SupLTS
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # SpinyStellate
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55], # TuftedIB
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55], # TuftedRS
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # DeepBasket
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # DeepAxoaxonic
                [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,26,27,31,32,33,34,35,36,37,38,39,40,44,45,46,47,48,49,50,51,52,53], # DeepLTS
                [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,38,39,40,41,42,43,44], #NontuftedRS
                [],
                []],
            [ # NontuftedRS
                [41,42,43,44], # SupPyrRS
                [41,42,43,44], # SupPyrFRB
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SupLTS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # SpinyStellate
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedIB
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47], # TuftedRS
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepBasket
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepAxoaxonic
                [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,44,45,46,47,48,49], # DeepLTS
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],  # NontuftedRS
                [6,7,8,9,10,11,12,13,14,19,20,21,22,23,24,25,26,27,32,33,34,35,36,37,38,39,40,45,46,47,48,49,50,51,52,53,58,59,60,61,62,63,64,65,66,71,72,73,74,75,76,77,78,79,84,85,86,87,88,89,90,91,92,97,98,99,100,101,102,103,104,105,110,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131], # TCR
                [2,3,4,15,16,17,28,29,30,41,42,43]], # nRT
            [ # TCR
                [45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68], # SupPyrRS
                [45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68], # SupPyrFRB
                [2,3,4,15,16,17,28,29,30,41,42,43], # SupBasket
                [2,3,4,15,16,17,28,29,30,41,42,43], # SupAxoaxonic
                [], # SupLTS
                [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53], # SpinyStellate
                [47,48,49,50,51,52,53,54,55], # TuftedIB
                [47,48,49,50,51,52,53,54,55], # TuftedRS
                [2,3,4,15,16,17,28,29,30,41,42,43], # DeepBasket
                [2,3,4,15,16,17,28,29,30,41,42,43], # DeepAxoaxonic
                [],
                [40,41,42,43,44], # NontuftedRS
                [],
                [2,3,4,15,16,17,28,29,30,41,42,43]], # nRT
            [ # nRT
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [1,2,15,28,41,54,67,80,93,106,119], # TCR
                [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53]] # nRT
            ]
        
    def check_pre_post_ratio(self):
        """Check the pre-post ratio for each celltype pair"""
        with open('connmatrix.txt') as connmatrix_file:
            index = -1
            for line in connmatrix_file.readlines():
                index += 1
                if index == 0:
                    print line
                    continue
                line = line.strip()
                if line == '':
                    continue
                values = line.split(',')
                # print values
                for ii in range(14):
                    if self.pre_post_ratio[index-1][ii] != int(values[ii]):
                        config.LOGGER.debug('pre-post ratio mismatch for %s-%s' % (self.celltype[index-1], self.celltype[ii]))
                        return False
        return True
                    
    def check_allowed_comps(self):
        """Compare the allowed_comps list to the original allowedcomp.py"""
        if len(self.allowed_comps) != 14:
            config.LOGGER.debug('Allowed_comp list must be of length 14. Got: %d' % (len(self.allowed_comps)))
            return False
        
        for index in range(len(self.allowed_comps)):
            if len(self.allowed_comps[index]) != 14:
                config.LOGGER.debug('Each row must lists for 14 cell types. But has %d for %s' % (len(self.allowed_comps[index]), self.celltype[index]))
                return False
        for ii in range(14):
            if self.celltype[ii].endswith('Axoaxonic'):
                continue
            for jj in range(14):
                try:
                    orig = allowedcomp.ALLOWED_COMP[self.celltype[ii]][self.celltype[jj]]
                    for kk in range(len(self.allowed_comps[ii][jj])):
                        if self.allowed_comps[ii][jj][kk] != orig[kk]:
                            config.LOGGER.debug('Valued don\'t match for %s -> %s' % (self.celltype[ii], self.celltype[jj]))
                            return False
                except KeyError:
                    if  self.allowed_comps[ii][jj] != []:
                        config.LOGGER.debug('Expected an empty list. But got something else for %s -> %s' % (self.celltype[ii], self.celltype[jj]))
                        return False
        return True


    def check_tau_ampa(self):
        """Compare the tau_AMPA with older version with dict."""        
        for pre_index in range(len(self.celltype)):
            for post_index in range(len(self.celltype)):
                left = self.tau_ampa[pre_index][post_index]
                try:
                    right = synapse.TAU_AMPA[self.celltype[pre_index]][self.celltype[post_index]]
                except KeyError:
                    if left!= 0.0:
                        config.LOGGER.debug('Key not present in synapse.py:TAU_AMPA - %s, %s' % (self.celltype[pre_index], self.celltype[post_index]))
                        return False
                    else:
                        right = 0.0
                if not numpy.allclose([left], [right]):
                    config.LOGGER.debug('Values not equal: TAU_AMPA - %s, %s: %g <> %g' % (self.celltype[pre_index], self.celltype[post_index], left, right))
                    return False
        return True

                        
    def check_tau_nmda(self):
        """Compare the tau_nmda with older version with dict."""        
        for pre_index in range(len(self.celltype)):
            for post_index in range(len(self.celltype)):
                left = self.tau_nmda[pre_index][post_index]
                try:
                    right = synapse.TAU_NMDA[self.celltype[pre_index]][self.celltype[post_index]]
                except KeyError:
                    if left != 0.0:
                        config.LOGGER.debug('Key not present in synapse.py:TAU_NMDA - %s, %s' % (self.celltype[pre_index], self.celltype[post_index]))
                        return False
                    else:
                        right = 0.0
                if not numpy.allclose([left], [right]):
                    config.LOGGER.debug('Values not equal: TAU_NMDA - %s, %s: %g <> %g' % (self.celltype[pre_index], self.celltype[post_index], left, right))
                    return False
        return True
                    
    def check_tau_gaba(self):
        """Compare the tau_gaba with older version with dict."""        
        for pre_index in range(len(self.celltype)):
            for post_index in range(len(self.celltype)):
                left = self.tau_gaba[pre_index][post_index]
                try:
                    if self.celltype[pre_index] == 'nRT':
                        right = synapse.TAU_GABA_FAST[self.celltype[pre_index]][self.celltype[post_index]]
                    else:
                        right = synapse.TAU_GABA[self.celltype[pre_index]][self.celltype[post_index]]
                except KeyError:
                    if left != 0.0:
                        config.LOGGER.debug('Key not present in synapse.py:TAU_GABA - %s, %s' % (self.celltype[pre_index], self.celltype[post_index]))
                        return False
                    else:
                        right = 0.0
                if not numpy.allclose([left], [right]):
                    config.LOGGER.debug('Values not equal: TAU_GABA - %s, %s: %g <> %g' % (self.celltype[pre_index], self.celltype[post_index], left, right))
                    return False
        return True
        
        

class TraubNet(object):
    def __init__(self, celltype_graph_file, cell_graph_file=None, format='gml', scale=1.0):
        """scale -- specifies the scale by which the cell count and
        edgecount should be multiplied when generating cell-cell
        graph.
        
        """
        self.celltype_graph_file = celltype_graph_file
        self.cell_graph_file = cell_graph_file
        self.format = format
        self.scale = scale
        self.nRT_g_gaba_high = 2.1e-9
        self.nRT_g_gaba_low = 0.7e-9
        self.__celltype_graph = None
        self.__cell_graph = None

    def load_cell_graph(self):
        self.__cell_graph = ig.read(self.cell_graph_file, format=self.format)

    def setup(self):
        start = datetime.now()
        if self.celltype_graph_file is None:
            self.generate_celltype_graph()
        else:
            self.__celltype_graph = ig.read(self.celltype_graph_file, format=self.format)
        print self.__celltype_graph.summary()
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('celltype_graph read from file %s in %g' % 
                                     (celltype_graph_file, delta.seconds + 1e-6 * delta.microseconds))

        # If ``cell_graph_file`` is specified, then try reading it
        # directly, otherwise generate the cell-cell graph using the
        # celltype-graph.  

        if self.cell_graph_file is not None:
            self.scale = 1.0
            start = datetime.now()
            self.__cell_graph = ig.read(cell_graph_file, format=format)
            end = datetime.now()
            delta = end - start
            config.BENCHMARK_LOGGER.info('cell_graph read from file %s in %g' % (cell_graph_file, delta.seconds + 1e-6 * delta.microseconds))
        else:
            self.generate_cell_graph()
        # Default rescaling - as described in the paper

        print self.__cell_graph.summary()


    def generate_celltype_graph(self):
        """Generate the celltype connectivity graph from hardcoded TraubFullNetData object."""
        tn = TraubFullNetData()
        graph = ig.Graph(0, directed=True)
        graph.add_vertices(len(tn.celltype))
        edge_count = 0
        start_index = 0
        for celltype in graph.vs:
            celltype['label'] = tn.celltype[celltype.index]
            celltype['count'] = tn.cellcount[celltype.index]
            celltype['start_index'] = start_index
            start_index += celltype['count']
            for posttype in graph.vs:
                pre_post_ratio = tn.pre_post_ratio[celltype.index][posttype.index]
                if pre_post_ratio > 0:
                    graph.add_edges((celltype.index, posttype.index))
                    edge_count += 1
                    graph.es[edge_count-1]['weight'] = pre_post_ratio
                    graph.es[edge_count-1]['g_ampa'] = tn.g_ampa_baseline[celltype.index][posttype.index]
                    graph.es[edge_count-1]['g_nmda'] = tn.g_nmda_baseline[celltype.index][posttype.index]
                    graph.es[edge_count-1]['g_gaba'] = tn.g_gaba_baseline[celltype.index][posttype.index]
                    graph.es[edge_count-1]['tau_ampa'] = tn.tau_ampa[celltype.index][posttype.index]
                    graph.es[edge_count-1]['tau_nmda'] = tn.tau_nmda[celltype.index][posttype.index]
                    graph.es[edge_count-1]['tau_gaba'] = tn.tau_gaba[celltype.index][posttype.index]
                    graph.es[edge_count-1]['ps_comps'] = str(tn.allowed_comps[celltype.index][posttype.index])
                    if celltype['label'] == 'nRT':
                        if posttype['label'] == 'TCR':
                            graph.es[edge_count-1]['tau_gaba_slow'] = tn.nRT_TCR_tau_gaba_slow
                        elif posttype['label'] == 'nRT':
                            graph.es[edge_count-1]['tau_gaba_slow'] = tn.nRT_nRT_tau_gaba_slow
        self.__celltype_graph = graph
        return graph
                            
            
            
    def generate_cell_graph(self):
        """generate the cell-cell graph from the celltype graph.

        """

        # * Idea: How can we specify the rules of generating cell-cell
        # graph in the celltype-graph itself?
        #       Right now, the code implicitly reflects my prior 
        # knowledge of the network from reading of the paper. But 
        # what about a general rule based network?

        start = datetime.now()
        
        # ``cell_graph`` is generated from the celltype_graph.  We
        # scale the number of cells of each type by ``scale`` and for
        # each celltype, set the ``start_index`` specifying the index
        # of the first vertex of this celltype in vertex sequence of
        # the cell-graph. ``start_index`` attribute is used when
        # adding the edges.

        cell_graph = ig.Graph(0, directed=True)
        for celltype in self.__celltype_graph.vs:
            count = int(celltype['count'] * self.scale)
            start_index = len(cell_graph.vs)
            celltype['start_index'] = start_index
            cell_graph.add_vertices(count)
            cell_graph.vs[start_index: start_index + count]['type_index'] = [celltype.index] * count
            cell_graph.vs[start_index: start_index + count]['label'] = ['%s_%d' % (celltype['label'], index) for index in range(count)]

        # Adding the edges is the tricky bit and this is the most
        # important part of the network definition.  

        # The connectivity in the model is described in terms of the
        # number of presynaptic cells per post synaptic cell. This
        # information is saved as the edge weights of the
        # ``celltype_graph``. We call this ``pre_post_ratio``. For
        # each edge in the ``celltype_graph``, we add a set of edges
        # from the vertex sequence for the presynaptic cells to each
        # post-synaptic cell. For each post-synaptic cell, we select
        # randomly ``pre_post_ratio`` no. of indices which fall within
        # the range of vertex indices for the presynaptic cell
        # (``start_index`` to ``start_index + count``). Moreover, for
        # each such edge, we select a random compartment no. from the
        # list of allowed compartments on the post-synaptic cell.
        edge_count = 0
        for edge in self.__celltype_graph.es:
            pre = edge.source
            post = edge.target
            pre_celltype = self.__celltype_graph.vs[pre]
            post_celltype = self.__celltype_graph.vs[post]
            pre_start_index = pre_celltype['start_index']
            post_start_index = post_celltype['start_index']
            pre_count = int(pre_celltype['count'] * self.scale)
            post_count = int(post_celltype['count'] * self.scale)
            pre_post_ratio = edge['weight']
            ps_comps = numpy.array(eval(edge['ps_comps']), dtype=int)
            config.LOGGER.debug('Connecting populations: pre - %s(%d), post - %s(%d), pre_post_ratio - %d' % (pre_celltype['label'], pre_count, post_celltype['label'], post_count, pre_post_ratio))
            if (pre_post_ratio == 0) or (len(ps_comps) == 0):
                continue
            pre_cell_indices = numpy.random.randint(low=pre_start_index, 
                                                    high=(pre_start_index + pre_count), 
                                                    size=(post_count, pre_post_ratio))
            post_comp_indices = numpy.random.randint(low=0, 
                                                     high=len(ps_comps), 
                                                     size=(post_count, pre_post_ratio))

            edge_list = [(int(pre_cell_index), post_cell_index + post_start_index) 
                         for post_cell_index in range(post_count)
                         for pre_cell_index in pre_cell_indices[post_cell_index]]            

            edge_start = len(cell_graph.es)
            new_edge_count = len(edge_list)
            edge_count += new_edge_count
            cell_graph.add_edges(edge_list)
            new_edges = cell_graph.es.select(range(edge_start, edge_start+new_edge_count))
            new_edges['pretype'] = [pre_celltype['label']] * new_edge_count
            new_edges['posttype'] = [post_celltype['label']] * new_edge_count
            new_edges['ps_comp'] = post_comp_indices.flatten()
            new_edges['tau_ampa'] = [edge['tau_ampa'] ] * new_edge_count
            new_edges['g_ampa'] = [edge['g_ampa']] * new_edge_count
            new_edges['tau_nmda'] = [edge['tau_nmda']] * new_edge_count
            new_edges['g_nmda'] = [edge['g_nmda']] * new_edge_count
            new_edges['tau_gaba'] = [edge['tau_gaba']] * new_edge_count
            # nRT->TCR and nRT->nRT GABA-ergic synapses have a slow component
            if pre_celltype['label'] == 'nRT' and ( post_celltype['label'] == 'TCR' or post_celltype['label'] == 'nRT' ):
                new_edges['tau_gaba_slow'] = [edge['tau_gaba_slow']] * new_edge_count
            # nRT->TCR GABA-ergic synapses are a special case as the baseline conductance is uniformly distributed between 0.7 and 2.1 nS.
            if  pre_celltype['label'] == 'nRT' and post_celltype['label'] == 'TCR':
                g_gaba = numpy.random.random_sample(new_edge_count) * (self.nRT_g_gaba_high - self.nRT_g_gaba_low) + self.nRT_g_gaba_low
            else:
                g_gaba = [edge['g_gaba']] * new_edge_count
            new_edges['g_gaba'] = g_gaba
        print 'Edges:', edge_count
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Cellgraph generated in %g s' % (delta.seconds + 1e-6 * delta.microseconds))
        self.__cell_graph = cell_graph
        return cell_graph

    def save_celltype_graph(self, filename, format=None):
        ig.write(self.__celltype_graph, filename, format=format)
    
    def save_cell_graph(self, filename, format=None):
        ig.write(self.__cell_graph, filename, format)

    def rescale_syn_conductance(self, pretype, posttype, conductance_type, scale_factor):
        """Rescale the conductance for synaptic conductance of
        conductance_type between pretype and posttype cells.

        pretype -- label of presynaptic celltype - can be a list

        posttype -- label of postsynaptic celltype - can be a list same in length as pretype

        conductance_type -- key for the specific conductance between pre and post type, e.g. 'g_gaba' - can be a list of such keys.

        scale_factor -- factor by which the conductance should be multiplied

        """
        if isinstance(pretype, str):
            if isinstance(posttype, str):
                pretype = [pretype]
                posttype = [posttype]
            else:
                raise Exception('Both pretype and posttype must be either string or lists of same size.')
        if isinstance(conductance_type, str):
            conductance_type = [conductance_type]
        assert len(pretype) == len(posttype)
        for index in range(len(pretype)):
            for edge in self.__cell_graph.es.select(pretype=pretype[index]).select(posttype=posttype[index]):
                for conductance_name in conductance_type:
                    edge[conductance_name] *= scale_factor
            
    def setup_model(self, container):
        """
        setup the actual MOOSE model from the cell_graph.
        
        container -- container for the cells in the network.
        """
        if self.__cell_graph is None:
            raise Exception('cell_graph is empty. First call generate_celltype_graph and generate_cell_graph to instantiate it.')
        if isinstance(container, str):
            if not config.context.exists(container):
                container = moose.Neutral(container)
        else:
            if not isinstance(container, moose.PyMooseBase):
                raise Exception('Container must be a MOOSE object.')
        for vertex in self.__celltype_graph.vs:
            cell_class = eval(vertex['label'])
            for cell_vertex in self.__cell_graph.vs.select(type_index=vertex.index):
                cell = cell_class(cell_vertex['label'], container)
                config.LOGGER.debug('Created %s' % (cell.path))
           
           


def testTraubFullNetData():
    """
    Do some checks to verify backward compatibility.
    """
    tn = TraubFullNetData()
    print 'check pre_post_ratio:', tn.check_pre_post_ratio()
    print 'check_tau_ampa:', tn.check_tau_ampa()
    print 'check_tau_nmda:', tn.check_tau_nmda()
    print 'check_tau_gaba:', tn.check_tau_gaba()
    print 'check_allowed_comps:', tn.check_allowed_comps()

if __name__ == '__main__':
    ## commented out for testing TraubFullNetData
    scale = 1.0
    if len(sys.argv) > 1:
        scale = float(sys.argv[1])
    # network = TraubNet('nx_celltype_graph.gml', scale=scale)
    network = TraubNet(None, scale=scale)
    network.generate_celltype_graph()
    network.generate_cell_graph()
    network.setup_model('traubnet')
    network.save_celltype_graph('celltype_graph.gml', format='gml')
    network.save_cell_graph('cell_graph.gml', format='gml')
    #! commented out till here for testing TraubFullNetData !
#    testTraubFullNetData()

# 
# ig_traubnet.py ends here
