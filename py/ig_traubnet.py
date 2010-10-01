#!/usr/bin/env python
# ig_traubnet.py --- 
# 
# Filename: ig_traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Sep 16 16:19:39 2010 (+0530)
# Version: 
# Last-Updated: Fri Oct  1 17:23:46 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 313
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
# 
# 

# Code:

import sys
import os
from datetime import datetime
import numpy
import igraph as ig
import config

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
                          500,
                          100,
                          100,
                          100,
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

        self.tau_nmda = [[ 130.0e-3,    130e-3,  100e-3,  100e-3, 100e-3, 130e-3, 130e-3, 130e-3,   100e-3,   100e-3,   100e-3,    130e-3,    0.0,    0.0],
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
        self.nRT_tau_gaba_slow = [10e-3, 44.5e-3]

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

class TraubNet(object):
    def __init__(self, celltype_graph_file, cell_graph_file=None, format='gml', scale=1.0):
        """scale -- specifies the scale by which the cell count and
        edgecount should be multiplied when generating cell-cell
        graph.
        
        """
        start = datetime.now()
        self.__celltype_graph = ig.read(celltype_graph_file, format=format)
        print self.__celltype_graph.summary()
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('celltype_graph read from file %s in %g' % 
                                     (celltype_graph_file, delta.seconds + 1e-6 * delta.microseconds))
        self.scale = 1

        # If ``cell_graph_file`` is specified, then try reading it
        # directly, otherwise generate the cell-cell graph using the
        # celltype-graph.  

        if cell_graph_file is not None:
            start = datetime.now()
            self.__cell_graph = ig.read(cell_graph_file, format=format)
            end = datetime.now()
            delta = end - start
            config.BENCHMARK_LOGGER.info('cell_graph read from file %s in %g' % (cell_graph_file, delta.seconds + 1e-6 * delta.microseconds))
        else:
            self.scale = float(scale) # scale is used only if we are generating cell-cell graph from scratch
            self.__cell_graph = self._generate_cell_graph()

        print self.__cell_graph.summary()

    def _generate_cell_graph(self):
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
            new_edges['ps_comp'] = post_comp_indices.flatten()                
        print 'Edges:', edge_count
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Cellgraph generated in %g s' % (delta.seconds + 1e-6 * delta.microseconds))
        return cell_graph

    def save_celltype_graph(self, filename, format=None):
        ig.write(self.__celltype_graph, filename, format=format)
    
    def save_cell_graph(self, filename, format=None):
        ig.write(self.__cell_graph, filename, format)
        
if __name__ == '__main__':
    scale = 1.0
    if len(sys.argv) > 1:
        scale = float(sys.argv[1])
    network = TraubNet('nx_celltype_graph.gml', scale=scale)
    network.save_celltype_graph('celltype_graph.gml', format='gml')
    network.save_cell_graph('cell_graph.gml', format='gml')
# 
# ig_traubnet.py ends here
