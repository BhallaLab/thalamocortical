# traubnet.py --- 
# 
# Filename: traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Aug 10 15:45:05 2010 (+0530)
# Version: 
# Last-Updated: Tue Aug 10 22:24:04 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 125
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This replaces TraubNet class in connection.py. I am switching from
# strings and maps to networkx graphs and attributes of the nodes and
# the edges.
# 
# 

# Change log:
# 
# 2010-08-10 15:48:26 (+0530) - initiated development.
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
import matplotlib.pyplot as plt

import networkx as nx

import config
import synapse

class TraubNet(object):
    """
    Store network information in graphs.
    
    """
    def __init__(self, 
                 connmatrix_file='connmatrix.txt', 
                 allowedcomp_file='allowedcomp.txt', 
                 cellcount_file='cells.txt'):
        """
        connmatrix_file -- filename of the text file containing
        connectivity matrix. This should be a csv file with header
        containining the cell types followed by rows of
        integers. row[i][j] is the number of presynaptic cells of type
        header[i] that connect to each cell of type header[j].
        
        allowedcomp_file -- filename of csv file containing the list
        of allowed postsynaptic compartment numbers for each cell type
        pair. A row of this file should be of the form:
        presynaptic_celltype, postsynaptic_celltype, n1, n2, ...
        
        cellcount_file -- filename of a text file containing all the
        cell types along with the size of their population.

        """
        self.__celltype_graph = self._read_celltype_graph(connmatrix_file, format='csv', cellcount_file=cellcount_file)
        self.__cell_graph = None
        
    def _read_celltype_graph(self, 
                             connmatrix_file, 
                             format='gml', 
                             cellcount_file=None):
        """Load the celltype-to-celltype connectivity map from file
        and return a nested dictionary dict where dict[X][Y] is the
        number of presynaptic cells of type X connecting to each
        postsynaptic cell of type Y.

        connmatrix_file -- the path of the file containing
        the connectivity matrix.
        
        format -- format of the file. Can be csv (comma separated
        values) or gml.

        The first row in the csv file should be the column headers - which
        are the cell types. The connection matrix itself is a square
        matrix with entry[i][j] specifying the number of presynaptic
        cells per postsynaptic cell, where the presynaptic cells are
        of type header[i] and the postsynaptic cell is of type
        header[j]

        """
        celltype_graph = None
        cellcount_dict = {}
        if cellcount_file is not None:
            with open(cellcount_file, 'r') as popfile:
                reader = csv.reader(popfile)
                for line in reader:
                    if len(line) > 0:
                        cellcount_dict[line[0]] = int(line[1])

        if format == 'csv':
            celltype_graph = nx.DiGraph()
            for key, value in cellcount_dict.items():
                print key, value
                celltype_graph.add_node(key, count=value)
            reader = csv.reader(file(connmatrix_file))
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
                    celltype_graph.add_edge(pre, post, weight=value)
                    col += 1
        elif format == 'gml':
            celltype_graph = nx.read_gml(connmatrix_file)
        celltype_graph.graph['doc'] = 'Celltype-based connectivity data. count of node *n* is the number of cells of type *n* that are present in the model. weight of edge (a, b) is the number of cells of type *a* that connect to each cell of type *b*.'
        return celltype_graph

    def plot_celltype_graph(self):
        """Display the celltype connectivity graph 

        """
        try:
            pos = nx.graphviz_layout(self.__celltype_graph)
        except Exception, e:
            print e
            pos = nx.spring_layout(self.__celltype_graph)
        node_size_list = [self.__celltype_graph.node[vertex]['count'] * 10 for vertex in self.__celltype_graph]
        edge_weights = [edata['weight'] for u, v, edata in self.__celltype_graph.edges(data=True)]
        nx.draw(self.__celltype_graph,
                pos,
                alpha=0.4,
                with_labels=True,
                node_size = node_size_list,
                edge_color=edge_weights,
                edge_cmap=plt.cm.jet,
                edge_vmin=min(edge_weights),
                edge_vmax=max(edge_weights))
        plt.show()

    def save_celltype_graph(self, filename='celltype_conn.gml', format='gml'):
        """
        Save the celltype-to-celltype connectivity information in a file.
        
        filename -- path of the file to be saved.

        format -- format to save in. Using GML as GraphML support is
        not complete in NetworkX.  

        """
        if format == 'gml':
            nx.write_gml(self.__celltype_graph, filename)
        else:
            raise Exception('Only format supported is gml. Received: %s' % (format))
        print 'Saved celltype connectivity graph in', filename

def test():
    net = TraubNet()
    net.plot_celltype_graph()
    net.save_celltype_graph()

if __name__ == '__main__':
    test()

# 
# traubnet.py ends here
