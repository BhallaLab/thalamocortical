# trbnet.py --- 
# 
# Filename: trbnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Oct 11 17:52:29 2010 (+0530)
# Version: 
# Last-Updated: Mon Nov 15 16:41:54 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 670
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
# 4. Testing: I tried to write some test code: but I realized midway
#    that this is rather circular. I am generating the graph using
#    manually entered data (TraubFullNetData). I am trying to validate
#    the celltype-graph against the celltype graph I generated
#    earlier. So if there is an error it will keep propagating. The only
#    reasonable test is manual check.
#
# 5. How to save the cell-graph?  After experimenting with
#    igraph/networkx and gml, graphml, pickle, formats, I tend towards
#    using simple adjacency matrices/edgelists with hdf5. Thus I can
#    have one single file for all the data (network info as well as
#    spikes).
#
#    More specifically, I'll save the final(post-scaling) g between
#    each pair of cells connected through a synapse.
#

# Change log:
#
# 2010-10-26 14:13:00 (+0530) finished implementation of methods
# _generate_celltype_graph and _read_celltype_graph
# 
# 2010-10-27 09:50:17 (+0530) implemented a simple test function to
# cross verify the celltype graph.
#
# 2010-11-09 17:44:03 (+0530) introduced pysparse module to utilize
# ll_mat for adjacency matrices for storing /g/.
# added creation of these matrices as part of generate_cell_graph fn.
# Updated scale_condunctance to take care of ggaba(nRT->TCR).

# Code:

from collections import defaultdict
from datetime import datetime
import igraph as ig
import numpy
import tables

from pysparse.sparse.spmatrix import ll_mat
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

class CellType(tables.IsDescription):
    """The CellType class is data for each row in the celltype table.
    name -- human readable name of the celltype

    index -- index of the celltype in the table. Do we need explicit index??

    count -- number of cells of this celltype.

    """
    name = tables.StringCol(16)
    index = tables.UInt8Col()
    count = tables.UInt16Col()

class TraubNet(object):
    """Implements the full network in Traub et al 2005 model.

    celltype_file -- A file containing the celltype-celltype-graph.

    graph_format -- format of the celltype and other graphs to be read or saved.

    scale -- scale factor for all the populations.        

    populations -- a dictionary mapping each cell-class-name to a list
    of cells of that class.

    index_cell_map -- dictionary mapping a global index to a cell instance

    cell_index_map -- dictionary mapping a cell instance to a global index

    """
    def __init__(self, celltype_file=None, format=None, scale=None):
        """
        celltype_file -- A file containing the celltype-celltype-graph.

        format -- format of the celltype and other graphs to be read or saved.

        scale -- scale factor for all the populations.
        """
        self.celltype_file = celltype_file
        self.scale = scale
        self.graph_format = format
        self.populations = defaultdict(list)
        self.celltype_graph = None
        self.cell_graph = None
        self.g_gaba_mat = None
        self.g_ampa_mat = None
        self.g_nmda_mat = None
        self.ps_comp_mat = None
        self.index_cell_map = {}
        self.cell_index_map = {}
        self.network_container = moose.Neutral('/net')
    
    def setup_from_celltype_file(self, celltype_file=None, format=None, scale=None):
        """Set up the network from a celltype-celltype graph file.

        celltype_file -- the file containing the celltype-celltype-graph.

        format -- format of the celltype_file

        scale -- scale factor for the network
        """
        if self.cell_graph is not None:
            del self.cell_graph
            self.cell_graph = None
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
        self.nRT_TCR_ggaba_low = tn.nRT_g_gaba_low
        self.nRT_TCR_ggaba_high = tn.nRT_g_gaba_high
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

    def _generate_cell_graph(self):
        start = datetime.now()
        self.cell_graph = ig.Graph(0, directed=True)
        total_count = 0
        
        for celltype in self.celltype_graph.vs:
            celltype['startindex'] = total_count
            cell_count = int(celltype['count'])
            cell_class = eval(celltype['label'])
            for ii in range(cell_count):
                cell = cell_class(cell_class.prototype, '%s/%s_%d' % (self.network_container.path, celltype['label'], ii))
                self.index_cell_map[total_count + ii] = cell
                self.cell_index_map[cell.id] = total_count + ii
            total_count += cell_count

        self.g_gaba_mat = ll_mat(total_count, total_count)
        self.g_ampa_mat = ll_mat(total_count, total_count)
        self.g_nmda_mat = ll_mat(total_count, total_count)
        self.ps_comp_mat = ll_mat(total_count, total_count)
        for edge in self.celltype_graph.es:
            pre = edge.source
            post = edge.target
            pretype = self.celltype_graph.vs[pre]
            posttype = self.celltype_graph.vs[post]
            prestart = int(pretype['startindex'])
            poststart = int(posttype['startindex'])
            precount = int(pretype['count'])
            postcount = int(posttype['count'])
            connprob = float(edge['weight'])
            ps_comps = numpy.array(eval(edge['pscomps']))
            config.LOGGER.debug('Connecting populations: pre=%s[:%d], post=%s[:%d], probability=%g' % (pretype['label'], pretype['count'], posttype['label'], posttype['count'], connprob))
            if connprob <= 0 or len(ps_comps) == 0:
                continue
            # pre_indices[i] is the array of global indices of the
            # presynaptic cells connecting to the i-th postsynaptic
            # cell of posttype.
            pre_indices = numpy.random.randint(low=prestart, high=prestart+precount, size=(postcount, int(connprob * precount)))
            # comp_indices[i][j] is the index of the postsynaptic
            # compartment in ps_comps for i-th postsynaptic
            # compartment for j-th presynaptic cell connecting to
            # postsynaptic cell
            comp_indices = numpy.random.randint(low=0, high=len(ps_comps), size=(postcount, int(connprob * precount)))
            # syn_list is the list of global index pairs for synapses
            syn_list = numpy.array([[preindex, postindex + poststart]
                                    for postindex in range(postcount)
                                    for preindex in pre_indices[postindex]],
                                   dtype='int32')
            self.ps_comp_mat.put(ps_comps[comp_indices.flatten()], syn_list[:,0], syn_list[:, 1])
            self.g_ampa_mat.put(float(edge['gampa']),
                                syn_list[:, 0], syn_list[:,1])
            self.g_nmda_mat.put(float(edge['gnmda']),
                                syn_list[:, 0], syn_list[:,1])
            if pretype['label'] == 'nRT' and posttype['label'] == 'TCR':
                self.g_gaba_mat.put(numpy.random.random_sample(len(syn_list)) * (self.nRT_TCR_ggaba_high - self.nRT_TCR_ggaba_low) + self.nRT_TCR_ggaba_low,
                                    syn_list[:,0],
                                    syn_list[:,1])
            else:
                self.g_gaba_mat.put(float(edge['ggaba']),
                                    syn_list[:,0], syn_list[:,1])
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('cell-cell network generation in: %g s' % (delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))

    def scale_populations(self, scale):
        """Scale the number of cells in each population by a factor."""
        if self.cell_graph is not None:
            raise Warning('Cell-graph already instantiated. Cannot rescale.')
        for vertex in self.celltype_graph.vs:
            vertex['count'] *= scale
    
    def scale_conductance(self, conductance_name, scale_dict):
        """Scale a particular synaptic conductance between pairs of celltypes.

        conductance_name -- key for the conductance to be scaled, e.g., 'ggaba'

        scale_dict -- a dict mapping celltype-pairs (pre, post) to scale factor. example:
                      {('SupPyrRS', 'SpinyStellate'): 0.5,
                       ('SupPyrFRB', 'SupLTS'): 1.5}
                       If either member of the tuple is '*', all celltypes are assumed.
        """
        for celltype_pair, scale_factor in scale_dict.items():
            if not isinstance(celltype_pair, tuple): 
                raise Warning('The keys in the scale_dict must be tuples of celltypes. Got %s of type %s instead' % (str(celltype_pair), type(celltype_pair)))
            if celltype_pair[0] == '*':
                pretype_seq = self.celltype_graph.vs
            else:
                pretype_seq = self.celltype_graph.vs.select(label_eq=celltype_pair[0])
            if celltype_pair[1] == '*':
                posttype_seq = self.celltype_graph.vs
            else:
                posttype_seq = self.celltype_graph.vs.select(label_eq=celltype_pair[1])
            for pretype in pretype_seq:
                for posttype in posttype_seq:
                    try:
                        edge_id = self.celltype_graph.get_eid(pretype.index, posttype.index)
                        edge = self.celltype_graph.es[edge_id]
                        if conductance_name in edge.attribute_names():
                            if conductance_name == 'ggaba' and pretype['label'] == 'nRT' and posttype['label'] == 'TCR':
                                self.nRT_TCR_ggaba_low *= scale_factor
                                self.nRT_TCR_ggaba_high *= scale_factor
                                edge[conductance_name] = 'uniform %f %f' % (self.nRT_TCR_ggaba_low, self.nRT_TCR_ggaba_high)
                            else:
                                edge[conductance_name] *= scale_factor

                    except Exception, e:
                        pass

    def _instantiate_model(self):
        for pretype in self.celltype_graph.vs:
            for posttype in self.celltype_graph.vs:                
                pass
        raise Exception('TODO: finish implementing this')


def test_generate_celltype_graph(celltype_file='celltype_graph.gml', format='gml'):
    celltype_graph = ig.read(celltype_file, format=format)
    trbnet = TraubNet()
    trbnet._generate_celltype_graph()
    for vertex in trbnet.celltype_graph.vs:
        original_vertex = celltype_graph.vs.select(label_eq=vertex['label'])
        assert original_vertex 
        assert original_vertex[0]['count'] == vertex['count']

    for edge in trbnet.celltype_graph.es:
        source = trbnet.celltype_graph.vs[edge.source]
        target = trbnet.celltype_graph.vs[edge.target]
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
            
            
def test_scale_conductance():
    netdata = TraubFullNetData()
    trbnet = TraubNet()
    trbnet._generate_celltype_graph()
    scale_dict_ampa = {('*', '*'): 2.0}
    scale_dict_nmda = {('*', 'SupLTS'): 0.2,
                       ('*', 'SupBasket'): 0.2,
                       ('*', 'SupAxoaxonic'): 0.2,
                       ('*', 'DeepLTS'): 0.2,
                       ('*', 'DeepBasket'): 0.2,
                       ('*', 'DeepAxoaxonic'): 0.2,
                       ('*', 'TCR'): 0.2,
                       ('*', 'nRT'): 0.2}
    
    scale_dict_ampa_ss_low = {('SpinyStellate', 'SpinyStellate'): 0.25/2.0}
    trbnet.scale_conductance('gampa', scale_dict_ampa)
    trbnet.scale_conductance('gnmda', scale_dict_nmda)
    trbnet.scale_conductance('gampa', scale_dict_ampa_ss_low)
    # TODO: now compare g's for each edge with the expected value
    for edge in trbnet.celltype_graph.es:
        source = trbnet.celltype_graph.vs[edge.source]
        target = trbnet.celltype_graph.vs[edge.target]
        src_index = netdata.celltype.index(source['label'])
        tgt_index = netdata.celltype.index(target['label'])
        gampa_baseline = netdata.g_ampa_baseline[src_index][tgt_index]
        gnmda_baseline = netdata.g_nmda_baseline[src_index][tgt_index]
        low_nmda_posttypes = [x[1] for x in scale_dict_nmda.keys()]
        if source['label'] == 'SpinyStellate' and target['label'] == 'SpinyStellate':
            assert numpy.allclose([edge['gampa']], [gampa_baseline * 0.25])
        else:
            assert numpy.allclose([edge['gampa']], [gampa_baseline * 2.0])
        if target['label'] in low_nmda_posttypes:
            assert numpy.allclose([edge['gnmda']], [gnmda_baseline * 0.2])
        else:
            assert numpy.allclose([edge['gnmda']], [gnmda_baseline])
    print 'test_scale_conductance: Successfully tested.'
        

if __name__ == '__main__':
    net = TraubNet()
    net._generate_celltype_graph()
    net._generate_cell_graph()

# 
# trbnet.py ends here
