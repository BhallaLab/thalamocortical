# trbnet.py --- 
# 
# Filename: trbnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Oct 11 17:52:29 2010 (+0530)
# Version: 
# Last-Updated: Wed Oct 20 22:17:40 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 72
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This is for putting together the lessons learnt from experiments
# with Graph as data structure to represent the netwrok.
# 
# I note the following points: 
#
# 1. It is nice to have the defining data-structure as a graph.  Each
#    node represents a homogeneous population. The node attributes can
#    store population specific information.
#
#    The connectivity information can be stored in the edges. The edge
#    attributes represent the information on the synapses. Edge weight
#    can represent pre-post ratio (no. of presynaptic cell of type
#    source vertex connecting to each cell of type destination vertex).
#
# 2. For analyzing the simulation data, it will be useful to have
# connectivity information for all relevant cells.

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

import config
import moose
from collections import defaultdict
import igraph as ig


class Traubnet(object):
    """Implements the full network in Traub et al 2005 model.

    cell_class -- name of the class of cells contained in this population.
    
    cells -- list containing the MOOSE cell objects.

    conn_graphs -- graphs for 
    """
    def __init__(self, celltype_file=None, scale=None):
        self.celltype_graph_file = celltype_file
        self.cells = None
        self.scale = scale
        self.populations = 

    def setup_populations(self):
        
    def setup(self):
        """TODO -- actually setup the connectivity here - maintain a
        cell-cell connection graph.
        """
        
        
    


# 
# trbnet.py ends here
