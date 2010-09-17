#!/usr/bin/env python
# ig_traubnet.py --- 
# 
# Filename: ig_traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Sep 16 16:19:39 2010 (+0530)
# Version: 
# Last-Updated: Fri Sep 17 17:04:08 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 221
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# igraph version of cell-cell graph generation for traub 2005 model
# 
# 

# Change log:
# 
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

class TraubNet(object):
    def __init__(self, celltype_graph_file, cell_graph_file=None, format='gml', scale=1.0):
        """scale -- specifies the scale by which the cell count and
        edgecount should be multiplied when generating cell-cell
        graph.
        
        """
        start = datetime.now()
        self.__celltype_graph = ig.read(celltype_graph_file, format=format)
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
            # config.LOGGER.debug('pre cell indices: low = %d, high = %s, size = (%d, %d)' % 
            #                     (pre_start_index, 
            #                      pre_start_index + pre_count, 
            #                      post_count, 
            #                      pre_post_ratio))
            post_comp_indices = numpy.random.randint(low=0, 
                                                     high=len(ps_comps), 
                                                     size=(post_count, pre_post_ratio))
            edge_list = [(int(pre_cell_index), post_cell_index + post_start_index) 
                         for post_cell_index in range(post_count)
                         for pre_cell_index in pre_cell_indices[post_cell_index]]
            edge_start = len(cell_graph.es)
            new_edge_count = len(edge_list)
            cell_graph.add_edges(edge_list)
            new_edges = cell_graph.es.select(range(edge_start, edge_start+new_edge_count))
            new_edges['ps_comp'] = post_comp_indices.flatten()                
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Cellgraph generated in %g s' % (delta.seconds + 1e-6 * delta.microseconds))
        return cell_graph
        
if __name__ == '__main__':
    scale = 1.0
    if len(sys.argv) > 1:
        scale = float(sys.argv[1])
    network = TraubNet('nx_celltype_graph.gml', scale=scale)
# 
# ig_traubnet.py ends here
