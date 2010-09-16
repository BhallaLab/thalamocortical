#!/usr/bin/env python
# ig_traubnet.py --- 
# 
# Filename: ig_traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Sep 16 16:19:39 2010 (+0530)
# Version: 
# Last-Updated: Fri Sep 17 00:28:37 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 151
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
        config.BENCHMARK_LOGGER.info('celltype_graph read from file %s in %g' % (celltype_graph_file, delta.seconds + 1e-6 * delta.microseconds))
        self.scale = 1
        if cell_graph_file is not None:
            start = datetime.now()
            self.__cell_graph = ig.read(cell_graph_file, format=format)
            end = datetime.now()
            delta = end - start
            config.BENCHMARK_LOGGER.info('cell_graph read from file %s in %g' % (cell_graph_file, delta.seconds + 1e-6 * delta.microseconds))
        else:
            self.scale = float(scale)
            self.__cell_graph = self._generate_cell_graph()

    def _generate_cell_graph(self):
        """generate the cell-cell graph from the celltype graph.
        """
        start = datetime.now()
        self.__cell_graph = ig.Graph(0, directed=True)
        for celltype in self.__celltype_graph.vs:
            count = int(celltype['count'] * self.scale)
            start_index = len(self.__cell_graph.vs)
            celltype['start_index'] = start_index
            self.__cell_graph.add_vertices(count)
            self.__cell_graph.vs[start_index: start_index + count]['type_index'] = [celltype.index] * count
            
        for edge in self.__celltype_graph.es:
            pre = edge.source
            post = edge.target
            pre_celltype = self.__celltype_graph.vs[pre]
            post_celltype = self.__celltype_graph.vs[post]
            pre_start_index = pre_celltype['start_index']
            post_start_index = post_celltype['start_index']
            pre_count = int(pre_celltype['count'] * self.scale)
            post_count = int(post_celltype['count'] * self.scale)
            print post_count
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
            edge_list = [(pre_cell_index, post_cell_index + post_start_index) 
                         for post_cell_index in range(post_count)
                         for pre_cell_index in pre_cell_indices[post_cell_index]]
            print 'Edges:', len(edge_list)
            for entry in edge_list:
                if len(entry) == 2 and isinstance(entry[0], int) and isinstance(entry[1], int):
                    print entry
                else:
                    print 'Not integer pairs:', entry
            edge_start = len(self.__cell_graph.es)
            new_edge_count = len(edge_list)
            self.__cell_graph.add_edges(edge_list)
            new_edges = self.__cell_graph.es.select(range(edge_start, edge_start+new_edge_count))
            new_edges['ps_comp'] = post_comp_indices.flat()                
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Cellgraph generated in %g s' % (delta.seconds + 1e-6 * delta.microseconds))
        
if __name__ == '__main__':
    scale = 1.0
    if len(sys.argv) > 1:
        scale = float(sys.argv[1])
    cell_graph = TraubNet('nx_celltype_graph.gml', scale=scale)
# 
# ig_traubnet.py ends here
