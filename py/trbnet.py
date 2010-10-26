# trbnet.py --- 
# 
# Filename: trbnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Oct 11 17:52:29 2010 (+0530)
# Version: 
# Last-Updated: Wed Oct 27 01:58:04 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 373
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
#    connectivity information for all relevant cells.
#
# 3. Decided to use HDF5 for data.  Need a clean way to associate
#    results of simulation with the model instance that was used. For
#    each cell's spike data, it should be possible to go back to the
#    model and see what was the connectivity for this cell. Thus model
#    is part of the data.
#
# 4. I tried to write some test code: but I realized midway that this
#    is rather circular. I am generating the graph using manually
#    entered data (TraubFullNetData). I am trying to validate the
#    celltype-graph against the celltype graph I generated earlier. So
#    if there is an error it will keep propagating. The only
#    reasonable test is manual check.
#

# Change log:
#
# 2010-10-26 14:13:00 (+0530) finished implementation of methods
# _generate_celltype_graph and _read_celltype_graph
# 
# 
# 

# Code:

from collections import defaultdict
import igraph as ig
import config
import moose

from trbnetdata import TraubFullNetData

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


# Is it a good idea to have a separate population class? What is the
# use?  When it comes to connectivity, the connectivity information is
# one level above the population. If a population is self contained,
# it should not know about other populations and thus not have
# connection probabilities dependent on the celltype of the other
# population.  If that connection information is not part of the
# population, it becomes just another container class. A set will do
# as well.

# class Population(object):
#     """Class to implement a homogeneous population"""


class TraubNet(object):
    """Implements the full network in Traub et al 2005 model.

    celltype_file -- A file containing the celltype-celltype-graph.

    graph_format -- format of the celltype and other graphs to be read or saved.

    scale -- scale factor for all the populations.        

    populations -- a dictionary mapping each cell-class-name to a list
    of cells of that class.

    """
    def __init__(self, celltype_file=None, format=None, scale=None):
        self.celltype_file = celltype_file
        self.scale = scale
        self.graph_format = format
        self.populations = defaultdict(list)
        self.celltype_graph = None
        self.cellgraph = None
    
    def setup_from_celltype_file(self, celltype_file=None, format=None, scale=None):
        """Set up the network from a celltype-celltype graph file.

        celltype_file -- the file containing the celltype-celltype-graph.

        format -- format of the celltype_file

        scale -- scale factor for teh network
        """
        if self.cellgraph is not None:
            return
        if celltype_file is not None:
            self.celltype_file = celltype_file
        if scale is not None:
            self.scale = scale
        elif self.scale is None:
            self.scale = 1.0
        if format is not None:
            self.format = format
        self._generate_celltype_graph()
                
    def setup(self):
        """Set up the master graph.
        """
        if self.celltype_file is None:
            config.LOGGER.info('Setting up network from predefined structure with full network information.')
            self._generate_celltype_graph()
        else:
            self._read_celltype_graph()
        
    def _generate_celltype_graph(self):
        """Generate the celltype-graph as in Traub model (without the
        gapjunctions at this point).

        The generated graph will be manipulated to control the model
        to be instantiated.

        """
        tn = TraubFullNetData()
        self.celltype_graph = ig.Graph(0, directed=True)
        self.celltype_graph.add_vertices(len(tn.celltype))
        edge_count = 0
        start_index = 0
        for celltype in self.celltype_graph.vs:
            celltype['label'] = tn.celltype[celltype.index]
            celltype['count'] = tn.cellcount[celltype.index]
            celltype['start_index'] = start_index
            start_index += celltype['count']
            for posttype in self.celltype_graph.vs:
                pre_post_ratio = tn.pre_post_ratio[celltype.index][posttype.index]
                if pre_post_ratio > 0:
                    self.celltype_graph.add_edges((celltype.index, posttype.index))
                    self.celltype_graph.es[edge_count]['weight'] = 1.0 * pre_post_ratio / celltype['count']
                    self.celltype_graph.es[edge_count]['gampa'] = tn.g_ampa_baseline[celltype.index][posttype.index]
                    self.celltype_graph.es[edge_count]['gnmda'] = tn.g_nmda_baseline[celltype.index][posttype.index]
                    self.celltype_graph.es[edge_count]['tauampa'] = tn.tau_ampa[celltype.index][posttype.index]
                    self.celltype_graph.es[edge_count]['taunmda'] = tn.tau_nmda[celltype.index][posttype.index]
                    self.celltype_graph.es[edge_count]['taugaba'] = tn.tau_gaba[celltype.index][posttype.index]
                    self.celltype_graph.es[edge_count]['pscomps'] = str(tn.allowed_comps[celltype.index][posttype.index])
                    self.celltype_graph.es[edge_count]['ekgaba'] = tn.ek_gaba[posttype.index]
                    if celltype['label'] == 'nRT':
                        if posttype['label'] == 'TCR':
                            self.celltype_graph.es[edge_count]['taugabaslow'] = tn.nRT_TCR_tau_gaba_slow
                        elif posttype['label'] == 'nRT':
                            self.celltype_graph.es[edge_count]['taugabaslow'] = tn.nRT_nRT_tau_gaba_slow
                            self.celltype_graph.es[edge_count]['ggaba'] = 'uniform %f %f' % (tn.nRT_g_gaba_low, tn.nRT_g_gaba_high) # How to specify distribution?
                    else:
                        self.celltype_graph.es[edge_count]['ggaba'] = tn.g_gaba_baseline[celltype.index][posttype.index]
                    edge_count += 1
                    
    def _read_celltype_graph(self):
        """
        read the celltype graph from a graph file.
        """
        self.celltype_graph = ig.read(self.celltype_file, format=self.graph_format)

    def test_generate_celltype_graph(self, celltype_file='celltype_graph.gml', format='gml'):
        celltype_graph = ig.read(celltype_file, format=format)
        self._generate_celltype_graph()
        for vertex in self.celltype_graph.vs:
            original_vertex = celltype_graph.vs.select(label_eq=vertex['label'])
            assert original_vertex 
            assert original_vertex[0]['count'] == vertex['count']

        for edge in self.celltype_graph.es:
            source = self.celltype_graph.vs[edge.source]
            target = self.celltype_graph.vs[edge.target]
            original_source = celltype_graph.vs.select(label_eq=source['label'])
            original_target = celltype_graph.vs.select(label_eq=target['label'])
            original_edge_id = celltype_graph.get_eid(original_source[0].index, original_target[0].index)
            original_edge = celltype_graph.es[original_edge_id]
            assert int(edge['weight'] * source['count']) == original_edge['weight']
            assert edge['gampa'] == original_edge['gampa']
            assert edge['gnmda'] == original_edge['gnmda']
            assert edge['tauampa'] == original_edge['tauampa']
            assert edge['taunmda'] == original_edge['taunmda']
            assert edge['taugaba'] == original_edge['taugaba']
            assert edge['pscomps'] == original_edge['pscomps']
            assert edge['ekgaba'] == original_edge['ekgaba']
            if source['label'] == 'nRT' and target['label'] == 'TCR':
                assert edge['taugabaslow'] == original_edge['taugabaslow']
            
            
                
        

if __name__ == '__main__':
    network = TraubNet()
    # network.setup()
    network.test_generate_celltype_graph()
    
# 
# trbnet.py ends here
