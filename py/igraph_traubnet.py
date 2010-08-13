# igraph_traubnet.py --- 
# 
# Filename: igraph_traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Aug 12 23:29:34 2010 (+0530)
# Version: 
# Last-Updated: Sat Aug 14 00:18:32 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 217
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This code is more or less same as traubnet.py - only uses igraph
# instead of networkx. The purpose is to compare the performance.
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

from datetime import datetime
import os
import csv
import numpy

import igraph
import config
import synapse

class TraubNet(object):
    def __init__(self,
                 connmatrix_file='connmatrix.txt', 
                 allowedcomp_file='allowedcomp.txt', 
                 cellcount_file='cells.txt'):
        self.__celltype_graph = self._read_celltype_graph(connmatrix_file, 
                                                          format='csv', 
                                                          cellcount_file=cellcount_file)
        self.__cell_graph = self._make_cell_graph('cell_graph.gml')
        print self.__cell_graph.summary()

    def _read_celltype_graph(self, 
                             connmatrix_file, 
                             format='gml', 
                             cellcount_file=None):
        celltype_graph = None
        cellcount_dict = {}
        if cellcount_file is not None:
            with open(cellcount_file, 'r') as popfile:
                reader = csv.reader(popfile)
                for line in reader:
                    if len(line) > 0:
                        cellcount_dict[line[0]] = int(line[1])

        if format == 'csv':
            celltype_graph = igraph.Graph(len(cellcount_dict), directed=True)
            index = 0
            for celltype, cellcount in cellcount_dict.items():
                celltype_graph.vs[index]['label'] = celltype 
                celltype_graph.vs[index]['count'] = cellcount
                index += 1
            reader = csv.reader(file(connmatrix_file))
            header = reader.next()
            row = 0
            for line in reader:
                if len(line) <= 0:
                    continue
                pre = header[row]
                source = celltype_graph.vs.select(label_eq=pre)
                if len(source) > 0:
                    source = source[0]
                row += 1
                col = 0
                for entry in line:
                    post = header[col]
                    target = celltype_graph.vs.select(label_eq=post)
                    if len(target) > 0:
                        target = target[0]
                    value = int(entry)
                    celltype_graph.add_edges((source.index, target.index))                    
                    eid = celltype_graph.get_eid(source.index, target.index)
                    celltype_graph.es[eid]['weight'] = value
                    col += 1
        elif format == 'gml':
            celltype_graph = igraph.read(f=connmatrix_file, format=format)

        celltype_graph['doc'] = 'Celltype-based connectivity data. \
count of node *n* is the number of cells of type *n* \
that are present in the model. weight of edge (a, b) \
is the number of cells of type *a* that connect to \
each cell of type *b*.'

        return celltype_graph

    def plot_celltype_graph(self):
        """Display the celltype connectivity graph 

        """
        return 

    def save_celltype_graph(self, filename='celltype_conn.gml', format='gml'):
        """
        Save the celltype-to-celltype connectivity information in a file.
        
        filename -- path of the file to be saved.

        format -- format to save in. Using GML as GraphML support is
        not complete in NetworkX.  

        """
        if format == 'gml':
            self.__celltype_graph.write(filename, format=format)
        else:
            raise Exception('Only format supported is gml. Received: %s' % (format))
        print 'Saved celltype connectivity graph in', filename

    def _make_cell_graph(self, filename=None):
        """Expand the celltype-to-celltype connectivity information
        and make a graph representing the network of the individual
        cells.

        Each cell is identified by the string {celltype}_{index}

        """
        if filename:
            try:
                start = datetime.now()
                cell_graph = igraph.read(filename, format='gml')
                end = datetime.now()
                delta = end - start
                config.BENCHMARK_LOGGER.info('Read cell_graph - time: %g s' 
                                             % (delta.seconds + 1e-6 * delta.microseconds))
                return cell_graph
            except Exception, e:
                print e, 'Now creating the cell_graph from scratch'
        start = datetime.now()
        cell_graph = igraph.Graph(0, directed=True)
        for celltype in self.__celltype_graph.vs:
            count = celltype['count']
            old_count = len(cell_graph.vs)
            cell_graph.add_vertices(count)
            for index in range(count):
                vertex = cell_graph.vs[index + old_count]
                vertex['type_index'] = celltype.index
                
        for edge in self.__celltype_graph.es:
            pre_type_index = edge.source
            post_type_index = edge.target
            post_cells = cell_graph.vs.select(type_index_eq=post_type_index)
            pre_cells = cell_graph.vs.select(type_index_eq=pre_type_index)
            pre_count = len(pre_cells)
            post_count = len(post_cells)
            pre_post_ratio = edge['weight']
            # randint returns unifrom random integers in [low, high)
            # interval. i-th row of pre_indices = list of indices of
            # presynaptic cells of type pre connected to i-th cell of
            # type post.
            pre_indices = numpy.random.randint(low=0, 
                                               high=pre_count, 
                                               size=(post_count, pre_post_ratio))
            edge_list = []
            for ii in range(post_count):
                post = post_cells[ii]    
                for jj in pre_indices[ii]:
                    pre = pre_cells[int(jj)]
                    edge_list.append((pre.index, post.index))
            cell_graph.add_edges(edge_list)

        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Built cell_graph programmatically - time: %g s' 
                                     % (delta.seconds + 1e-6 * delta.microseconds))
        return cell_graph

    def plot_cell_graph(self):
        """Display the cell-to-cell connection graph.

        """
        igraph.plot(self.__cell_graph)

    def save_cell_graph(self, filename='cell_graph.gml'):
        """Save the cell to cell connectivity graph in a GML file.

        """
        igraph.write(self.__cell_graph, filename, format='gml')
        print 'Saved cell-to-cell connectivity data in', filename
    
def test():
    net = TraubNet()
    # net.plot_celltype_graph()
    # net.save_celltype_graph()
    # net.plot_cell_graph()
    net.save_cell_graph()
    
if __name__ == '__main__':
    test()
                 


# 
# igraph_traubnet.py ends here
