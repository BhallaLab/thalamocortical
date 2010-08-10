# connection.py --- 
# 
# Filename: connection.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Jun 26 17:22:01 2010 (+0530)
# Version: 
# Last-Updated: Mon Aug  9 17:54:07 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 417
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This class is for setting up the connection
# information. It does not need the model itself to be instantiated,
# but expands the paths as strings and generates the random numbers to
# select the compartments.
# 
# 
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

from collections import defaultdict
from datetime import datetime
import os
import csv
import numpy
import networkx as nx
import igraph

import matplotlib.pyplot as plt

import synapse
import config


class TraubNet:
    """
    Store network information as strings and maps.

    cellcount -- dictionary with cell-type as key and cell-count as
    value.

    presynaptic -- dictionary of cell types and presynaptic
    compartment no.

    forward -- presynaptic cell no. mapped to set of synchans on
    postsynaptic compartment.

    backward -- postsynaptic synchans mapped to set of presynaptic
    cells.
    
    weight -- dictionary of dictionaries:
    weight[precell][post_synchan] = value of gbar for this connection.

    celltype -- list of cell types - lexicographically sorted

    cellnet -- a graph containing cell-to cell connectivity.

    compnet -- a graph containing compartment level connectivity.

    synnet -- a graph containing synchan level connectivity
    """
    def __init__(self, connmatrix_file='connmatrix.txt', allowedcomp_file='allowedcomp.txt', cellcount_file='cells.txt'):
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
        self.__connmap = None
        self.__allowedcompmap = None        
        self.__population = None
        self.__celltype = None
        self.__cellcount = {}
        self.cellcount_file = cellcount_file
        self.__connmap_file = connmatrix_file
        self.allowedcompmap_file = allowedcomp_file
        self.presynaptic = {}
        self.forward = defaultdict(set)
        self.backward = defaultdict(set)
        self.weight = defaultdict(dict)
        self.cellnet = nx.MultiDiGraph()
        self.cellclassnet = nx.DiGraph()
        self.i_cellclassnet = igraph.Graph(n=len(self.cellcount), directed=True) # igraph version of the cellclassnet
        self.cellclass_vertex_map = {}
        self.vertex_cellclass_map = {}
#        raise NotImplementedError('TODO: incorporate networkx module\'s MultiDiGraph class to implement cell to cell/comp/synchan connectivity.')
        
        self._expand()

    def _get_connmap(self):
        """Load the celltype-to-celltype connectivity map from file
        and return a nested dictionary dict where dict[X][Y] is the
        number of presynaptic cells of type X connecting to each
        postsynaptic cell of type Y.

        it can be set to filename -- the path of a csv file containing
        the connectivity matrix.

        The first row in the file should be the column headers - which
        are the cell types. The connection matrix itself is a square
        matrix with entry[i][j] specifying the number of presynaptic
        cells per postsynaptic cell, where the presynaptic cells are
        of type header[i] and the postsynaptic cell is of type
        header[j]

        """
        if self.__connmap is not None:
            return self.__connmap
        self.__connmap = defaultdict(dict)
        reader = csv.reader(file(self.__connmap_file))
        header = reader.next()
        row = 0
        for line in reader:
            if len(line) <= 0:
                continue
            pre = header[row]
            row = row + 1
            col = 0
            for entry in line:
                post = header[col]
                value = int(entry)
                self.__connmap[pre][post] = value
                col = col + 1
                
        return self.__connmap

    def _set_connmap(self, filename):
        if filename == self.__connmap_file:
            return
        if not os.path.exists(filename):
            raise Exception('%s does not exist.' % (filename))
        self.__connmap_file = filename
        self.__connmap = None

    # Define connmap as a property
    connmap = property(_get_connmap, _set_connmap)

    def _get_allowedcompmap(self):
        """Load the tables for allowed compartment list for synapses
        between each types of cells. 
        
        Return a nested dictionary with entries like:
        entry["X"]["Y"] = [n1, n2, n3, n4, ...]
        
        where only compartments with index n1, n2, n3, .... are
        allowed as postsynaptic compartment for synapses from cells of
        type X to cells of type Y.
        
        it can be set to filename -- the source of the map. This file
        contains source celltype, destination celltype followed by the
        allowed compartmnets.
        
        """
        if self.__allowedcompmap is not None:
            return self.__allowedcompmap
        with open(self.allowedcompmap_file, 'r') as cmpmap_file:
            try:
                reader = csv.reader(cmpmap_file)
                self.__allowedcompmap = defaultdict(dict)
                row = 0
                for line in reader:
                    if len(line) <= 0:
                        continue
                    self.__allowedcompmap[line[0]][line[1]] = []
                    for entry in line[2:]:
                        self.__allowedcompmap[line[0]][line[1]].append(int(entry))
            except csv.Error, e:
                print 'file %s, line %d: %s' % (filename, reader.line_num, e)
                raise 

        return self.__allowedcompmap

    def _set_allowedcompmap(self, filename):
        if filename == self.allowedcompmap_file:
            return
        if not os.path.exists(filename):
            raise Exception('%s does not exist.' % (filename))
        self.allowedcompmap_file = filename
        self.__allowedcompmap = None

    # Define allowedcompmap as a property
    allowedcompmap = property(_get_allowedcompmap, _set_allowedcompmap)

    @property
    def cellcount(self):
        if self.__cellcount:
            return self.__cellcount
        with open(self.cellcount_file, 'r') as popfile:
            reader = csv.reader(popfile)
            for line in reader:
                if len(line) > 0:
                    self.__cellcount[line[0]] = int(line[1])
        return self.__cellcount

    @property
    def celltype(self):
        """A list of cell types. It is sorted alphabetically."""
        if self.__celltype is not None:
            return self.__celltype
        self.__celltype = sorted(self.cellcount.keys())
        return self.__celltype

    def _expand(self):
        """
        expand the network into the maps.
        """
        config.LOGGER.debug(__name__ + ': starting.')
        index = 0
        min_count = min(self.cellcount.values())
        for cell in self.cellcount.keys():
            print '$$', cell
            self.cellclass_vertex_map[cell] = index
            self.vertex_cellclass_map[index] = cell
            self.i_cellclassnet.vs[index]['label'] = cell
            self.i_cellclassnet.vs[index]['size'] = self.cellcount[cell]*1.0/min_count
            index += 1
        for v in self.i_cellclassnet.vs:
            print '???', v
        print self.i_cellclassnet.vs
        self.min_prepost_count = 1000000
        self.max_prepost_count = 0
        for pretype in self.celltype:
            for posttype in self.celltype:
                precell_count = self.connmap[pretype][posttype]
                self.cellclassnet.add_node(pretype, count=self.cellcount[pretype])
                if precell_count > 0:
                    if precell_count > self.max_prepost_count:
                        self.max_prepost_count = precell_count
                    if precell_count < self.min_prepost_count:
                        self.min_prepost_count = precell_count
                    self.cellclassnet.add_edge(pretype, posttype, weight=precell_count, prepost_ratio=precell_count)
                    edge = (self.cellclass_vertex_map[pretype], self.cellclass_vertex_map[posttype])
                    self.i_cellclassnet.add_edges([edge])
                print '##', pretype, ':', precell_count, posttype
        # This is the unfriendly bit of igraph - edges have their own
        # existence, we cannot just set attribute on a pair of
        # vertices.
        for edge in self.i_cellclassnet.es:
            source_cell = self.vertex_cellclass_map[edge.source]
            target_cell = self.vertex_cellclass_map[edge.target]
            edge['prepost_ratio'] = self.connmap[source_cell][target_cell]
            # edge['arrow_size'] = self.connmap[source_cell][target_cell]*1.0/(10.0*self.min_prepost_count)
        
        config.LOGGER.debug(__name__ + ': finished.')

    def draw_cellclassnet(self):
        try:
            pos = nx.graphviz_layout(self.cellclassnet)
        except:
            pos = nx.graphviz_layout(self.cellclassnet)
        edge_widths = [edata['prepost_ratio']*1.0/self.min_prepost_count for u, v, edata in self.cellclassnet.edges(data=True)]
        print edge_widths
        node_sizes = [self.cellclassnet.node[v]['count']*10 for v in self.cellclassnet]
        
        # nx.draw_networkx_nodes(self.cellclassnet, pos, node_size=node_sizes, with_labels=True)
        # nx.draw_networkx_edges(self.cellclassnet, pos, 
        #                        alpha=0.4, 
        #                        width=edge_widths,
        #                        edge_color=edge_widths, 
        #                        edge_cmap=plt.cm.jet, 
        #                        edge_vmin=min(edge_widths), 
        #                        edge_vmax=max(edge_widths))

        nx.draw(self.cellclassnet, pos, 
                alpha=0.4,
                with_labels=True, 
                node_size=node_sizes,
                # width=5,
                # width=edge_widths,
                edge_color=edge_widths,
                edge_cmap=plt.cm.jet,
                edge_vmin=1.0,
                edge_vmax=1.0*self.max_prepost_count/self.min_prepost_count)
        # nx.draw_graphviz(self.cellclassnet, pos, alpha=0.4, with_labels=True, node_size=node_sizes, edge_color=edge_widths, edge_cmap=plt.cm.jet, edge_vmin=1.0, edge_vmax=1.0*self.max_prepost_count/self.min_prepost_count)
        plt.show()

    def i_draw_cellclassnet(self):
        print self.i_cellclassnet.vs['label']
        igraph.plot(self.i_cellclassnet, layout='fr', edge_color='blue')

    def gv_draw_cellclassnet(self):
        print 'Saving as a dot file.'
        # A_graph = nx.to_agraph(self.cellclassnet)
        nx.write_dot(self.cellclassnet, 'cellclassnet.dot')
        print 'Saved cellclassnet.dot'
                              

def test():
    net = TraubNet()
    net.i_draw_cellclassnet()
    net.draw_cellclassnet()
    net.gv_draw_cellclassnet()

if __name__ == '__main__':
    test()
# 
# connection.py ends here
