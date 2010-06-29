# connection.py --- 
# 
# Filename: connection.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Jun 26 17:22:01 2010 (+0530)
# Version: 
# Last-Updated: Tue Jun 29 10:06:41 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 138
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
import allowedcomp
import synapse
import config
class TraubNet:
    """
    Network information as strings and maps.

    population -- dictionary with cell-type as key and cell-count as
    value.

    presynaptic -- dictionary of cell types and presynaptic
    compartment no.

    forward -- presynaptic compartment mapped to set of synchans on
    the postsynaptic compartment.

    backward -- postsynaptic synchans mapped to set of presynaptic
    cells.
    
    weight -- dictionary of dictionaries:
    weight[precell][post_synchan] = value of gbar for this connection.
    """
    def __init__(self, connmatrix='connmatrix.txt', allowedcomp='allowedcomp.txt', cells='cells.txt'):
        """
        connmatrix -- filename of the text file containing
        connectivity matrix. This should be a csv file with header
        containining the cell types followed by rows of
        integers. row[i][j] is the number of presynaptic cells of type
        header[i] that connect to each cell of type header[j].

        allowedcomp -- filename of csv file containing the list of
        allowed postsynaptic compartment numbers for each cell type
        pair. A row of this file should be of the form:
        presynaptic_celltype, postsynaptic_celltype, n1, n2, ...

        """
        self.__connmap = None
        self.__allowedcompmap = None
        self.connmap_file = connmatrix
        self.allowedcompmap_file = allowedcomp
        self.population = {}
        self.presynaptic = {}
        self.forward = defaultdict(set)
        self.backward = defaultdict(set)
        self.weight = defaultdict(dict)
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
        reader = csv.reader(file(self.connmap_file))
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
        if filename == self.connmap_file:
            return
        if not os.path.exists(filename):
            raise Exception('%s does not exist.' % (filename))
        self.connmap_file = filename
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
        self.__allowedcompmap = defaultdict(dict)
        reader = csv.reader(file(self.allowedcompmap_file))
        row = 0
        for line in reader:
            if len(line) <= 0:
                continue
            self.__allowedcompmap[line[0]][line[1]] = []
            for entry in line[2:]:
                self.__allowedcompmap[line[0]][line[1]].append(int(entry))
        
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

    def _expand(self):
        """
        expand the network into the maps.
        """
        raise NotImplementedError
        
# 
# connection.py ends here
