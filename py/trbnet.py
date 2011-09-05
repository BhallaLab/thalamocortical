# trbnet.py --- 
# 
# Filename: trbnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Oct 11 17:52:29 2010 (+0530)
# Version: 
# Last-Updated: Mon Sep  5 10:34:18 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 1705
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
#
# 2010-12-28 14:26:47 (+0530) fixed serious mistake in get_maxdegree_cell_indices.


# Code:

from collections import defaultdict
from datetime import datetime
import igraph
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
import synapse


layer_to_celltype = {
    '1': [],
    '2/3': ['SupPyrRS', 'SupPyrFRB', 'SupBasket', 'SupAxoaxonic', 'SupLTS'],
    '4': ['SpinyStellate'],
    '5': ['TuftedIB', 'TuftedRS'],
    '6': ['NontuftedRS', 'DeepBasket', 'DeepAxoaxonic', 'DeepLTS'],
    'Thalamus': ['nRT', 'TCR']
}

celltype_to_layer = {
    'SupPyrRS': ['2/3'],
    'SupPyrFRB': ['2/3'],
    'SupBasket': ['2/3'],
    'SupAxoaxonic': ['2/3'],
    'SupLTS': ['2/3'],
    'SpinyStellate': ['4'],
    'TuftedIB': ['5'],
    'TuftedRS': ['5'],
    'NonTuftedRS': ['6'],
    'DeepBasket': ['6'],
    'DeepAxoaxonic': ['6'],
    'DeepLTS': ['6'],
    'nRT': ['Thalamus'],
    'TCR': ['Thalamus']
}
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

    This is the datastructure used for saving celltype info hdf5 format.
    
    name -- human readable name of the celltype

    index -- index of the celltype in the table. Do we need explicit index??

    count -- number of cells of this celltype.

    """
    name = tables.StringCol(16)
    index = tables.Int8Col()
    count = tables.Int16Col()

class SynEdge(tables.IsDescription):
    """Describes an edge of celltype graph. This is used for saving
    synase information in HDF5 file."""
    # I think I should just subclass instead of keeping multiple alternate data entries
    source = tables.Int8Col()    # Index of source celltype
    target =  tables.Int8Col()   # Index of destination celltype
    weight =  tables.Float64Col() # connection probability from sourec to target
    gampa = tables.Float64Col()   # max conductance for source -> target AMPA synapse
    gnmda =  tables.Float64Col()  # max conductance for source -> target NMDA synapse
    tauampa = tables.Float64Col() # decay time constant for source -> target AMPA synapse
    taunmda =  tables.Float64Col() # decay time constant for source -> target NMDA synapse
    tau2nmda =  tables.Float64Col() # rise time constant for source -> target NMDA synapse (from NEURON mod file)
    taugaba =  tables.Float64Col()  # decay time constant for source -> target GABA synapse
    taugabaslow =  tables.Float64Col() # slower decay time constant for source -> target GABA synapse
    pscomps =  tables.UInt8Col(shape=(90, ))    # compartment nos in target cell where synapses are allowed
    ekgaba = tables.Float64Col() # reversal potential for gaba synapses
    ggaba =  tables.Float64Col(shape=(2)) # gaba conductance (distributed uniformly between first and second entry)
    prelease = tables.Float64Col() # Baseline synaptic release probability for this pair

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
    def __init__(self, celltype_file=None, format=None, scale=None, container=None):
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
        self.electrodes = []
        self.electrode_container = moose.Neutral('/lfp')
        if container is None:
            self.network_container = moose.Neutral('/net')
        elif isinstance(container, str):
            self.network_container = moose.Neutral('%s/net', container)
        elif isinstance(container, moose.Neutral):
            self.network_container = moose.Neutral('net', container)
        else:
            raise('Need a moose-object/string/None as container. Got %s of type %s' % (container, container.__class__.__name__))
        self.ectopic_container = moose.Neutral('ectopic_spikes', self.network_container.parent)
        self.tweaks_doc = []
        
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
        self.celltype_graph = igraph.Graph(0, directed=True)
        self.celltype_graph.add_vertices(len(tn.celltype))
        self.nRT_TCR_ggaba_low = tn.nRT_TCR_ggaba_low
        self.nRT_TCR_ggaba_high = tn.nRT_TCR_ggaba_high
        self.frac_nRT_TCR_ggaba_fast = tn.frac_nRT_TCR_gaba_fast
        self.frac_nRT_nRT_ggaba_fast = tn.frac_nRT_nRT_gaba_fast
        edge_count = 0
        for celltype in self.celltype_graph.vs:
            celltype['label'] = tn.celltype[celltype.index]
            celltype['count'] = tn.cellcount[celltype.index]
            celltype['ectopicinterval'] = tn.ectopic_interval[celltype.index]
            for posttype in self.celltype_graph.vs:
                pre_post_ratio = tn.pre_post_ratio[celltype.index][posttype.index]
                if pre_post_ratio > 0:
                    self.celltype_graph.add_edges((celltype.index, posttype.index))
                    new_edge = self.celltype_graph.es[edge_count]
                    new_edge['weight'] = 1.0 * pre_post_ratio / celltype['count']
                    new_edge['gampa'] = tn.g_ampa_baseline[celltype.index][posttype.index] * tn.tau_ampa[celltype.index][posttype.index] / numpy.e # This is how gmax is related to c (baseline conductance scaling factor) for AMPA in traub model
                    new_edge['gnmda'] = tn.g_nmda_baseline[celltype.index][posttype.index]
                    new_edge['tauampa'] = tn.tau_ampa[celltype.index][posttype.index]
                    new_edge['taunmda'] = tn.tau_nmda[celltype.index][posttype.index]
                    new_edge['taugaba'] = tn.tau_gaba[celltype.index][posttype.index]
                    new_edge['pscomps'] = str(tn.allowed_comps[celltype.index][posttype.index])
                    new_edge['ekgaba'] = tn.ek_gaba[posttype.index]
                    new_edge['ggaba'] = tn.g_gaba_baseline[celltype.index][posttype.index]
                    new_edge['prelease'] = tn.p_release[celltype.index][posttype.index]
                    if celltype['label'] == 'nRT':
                        if posttype['label'] == 'TCR':
                            new_edge['taugabaslow'] = tn.nRT_TCR_tau_gaba_slow
                            new_edge['ggaba'] = 'uniform %g %g' % (tn.nRT_TCR_ggaba_low, tn.nRT_TCR_ggaba_high) # How to specify distribution?
                        elif posttype['label'] == 'nRT':
                            new_edge['taugabaslow'] = tn.nRT_nRT_tau_gaba_slow
                    
                    edge_count += 1

    def _read_celltype_graph(self):
        """
        read the celltype graph from a graph file.
        """
        self.celltype_graph = igraph.read(self.celltype_file, format=self.graph_format)

    def _generate_cell_graph(self):
        start = datetime.now()
        self.cell_graph = igraph.Graph(0, directed=True)
        total_count = 0
        for celltype in self.celltype_graph.vs:
            celltype['startindex'] = total_count
            cell_count = int(celltype['count'])
            config.LOGGER.debug('%s - population size %d' % (celltype['label'], cell_count))
            total_count += cell_count
        config.LOGGER.info('Total cell count: %d' % (total_count))
        self.g_gaba_mat = ll_mat(total_count, total_count)
        self.g_ampa_mat = ll_mat(total_count, total_count)
        self.g_nmda_mat = ll_mat(total_count, total_count)
        self.ps_comp_mat = ll_mat(total_count, total_count)
        for edge in self.celltype_graph.es:
            pre = edge.source
            post = edge.target
            pretype_vertex = self.celltype_graph.vs[pre]
            posttype_vertex = self.celltype_graph.vs[post]
            prestart = int(pretype_vertex['startindex'])
            poststart = int(posttype_vertex['startindex'])
            precount = int(pretype_vertex['count'])
            postcount = int(posttype_vertex['count'])
            connprob = float(edge['weight'])
            
            ps_comps = numpy.array(eval(edge['pscomps']), dtype=numpy.float)
            config.LOGGER.debug('Connecting populations: pre=%s[:%d], post=%s[:%d], probability=%g' % (pretype_vertex['label'], pretype_vertex['count'], posttype_vertex['label'], posttype_vertex['count'], connprob))
            config.LOGGER.debug('ggaba= %s, type:%s' % (str(edge['ggaba']), edge['ggaba'].__class__.__name__))
            config.LOGGER.debug('allowed postsynaptic compartments: %s (after conversion: %s)' % (edge['pscomps'], ps_comps))
            pre_per_post =  int(connprob * precount)
            if (connprob <= 0) or (len(ps_comps) == 0) or (precount <= 0) or (postcount <= 0) or (pre_per_post <= 0):
                continue
            # pre_indices[i] is the array of global indices of the
            # presynaptic cells connecting to the i-th postsynaptic
            # cell of posttype.
            pre_indices = numpy.random.randint(low=prestart, high=prestart+precount, size=(postcount,pre_per_post))
            # comp_indices[i][j] is the index of the postsynaptic
            # compartment in ps_comps for i-th postsynaptic
            # compartment for j-th presynaptic cell connecting to
            # postsynaptic cell
            comp_indices = numpy.random.randint(low=0, high=len(ps_comps), size=(postcount, pre_per_post))
            # syn_list is the list of global index pairs for synapses
            syn_list = numpy.array([[preindex, postindex + poststart]
                                    for postindex in range(postcount)
                                    for preindex in pre_indices[postindex]],
                                   dtype=numpy.int32)
            config.LOGGER.debug(edge['pscomps'])
            indices = comp_indices.flatten()
            ps_comp_list = ps_comps[indices]
            config.LOGGER.debug('ps_comps list has length: %d, syn_list has length: %d' % (len(ps_comp_list), len(syn_list)))
            self.ps_comp_mat.put(ps_comp_list, syn_list[:,0], syn_list[:, 1])
            self.g_ampa_mat.put(float(edge['gampa']),
                                syn_list[:, 0], syn_list[:,1])
            self.g_nmda_mat.put(float(edge['gnmda']),
                                syn_list[:, 0], syn_list[:,1])
            if (pretype_vertex['label'] == 'nRT') and (posttype_vertex['label'] == 'TCR'):
                self.g_gaba_mat.put(numpy.random.random_sample(len(syn_list)) * (self.nRT_TCR_ggaba_high - self.nRT_TCR_ggaba_low) + self.nRT_TCR_ggaba_low,
                                    syn_list[:,0],
                                    syn_list[:,1])
            else:
                self.g_gaba_mat.put(float(edge['ggaba']),
                                    syn_list[:,0], syn_list[:,1])                    

        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('cell-cell network generation in: %g s' % (delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))

    def create_network(self):
        """Instantiate the network in MOOSE"""
        config.LOGGER.debug('Creating network')
        synchan_classname = 'SynChan'
        nmdachan_classname = 'NMDAChan'
        if config.stochastic:
            synchan_classname = 'STPSynChan'
            nmdachan_classname = 'STPNMDAChan'
        starttime = datetime.now()
        total_count = 0
        for celltype in self.celltype_graph.vs:
            cell_count = int(celltype['count'])
            cell_class = eval(celltype['label'])
            for ii in range(cell_count):
                cell = cell_class(cell_class.prototype, '%s/%s_%d' % (self.network_container.path, celltype['label'], ii))
                self.index_cell_map[total_count + ii] = cell
                self.cell_index_map[cell] = total_count + ii
                self.populations[celltype['label']].append(total_count + ii)
            total_count += cell_count
        for syn_edge in self.celltype_graph.es:
            pretype_vertex = self.celltype_graph.vs[syn_edge.source]
            posttype_vertex = self.celltype_graph.vs[syn_edge.target]
            delay = synapse.SYNAPTIC_DELAY_DEFAULT
            p_release = syn_edge['prelease']
            if pretype_vertex['label'] in ['nRT', 'TCR']:            
                if posttype_vertex['label'] not in ['nRT', 'TCR']:
                    delay = synapse.SYNAPTIC_DELAY_THALAMOCORTICAL
            else:
                if posttype_vertex['label'] in ['nRT', 'TCR']:
                    delay = synapse.SYNAPTIC_DELAY_CORTICOTHALAMIC
            prestart = int(pretype_vertex['startindex'])
            poststart = int(posttype_vertex['startindex'])
            precount = int(pretype_vertex['count'])
            postcount = int(posttype_vertex['count'])
            for pre_index in range(prestart, prestart+precount):
                precell = self.index_cell_map[pre_index]
                precomp = precell.comp[precell.presyn]
                for post_index in range(poststart, poststart+postcount):
                    postcell = self.index_cell_map[post_index]
                    postcompindex = int(self.ps_comp_mat[pre_index, post_index])
                    if postcompindex > 255:
                        raise Exception('%s->%s -- PS comp has absurd index %d' % (precell.path, postcell.path, postcompindex))
                    postcomp = postcell.comp[postcompindex]
                    if postcomp is None:
                        continue
                    g_ampa = self.g_ampa_mat[pre_index, post_index]
                    if g_ampa != 0.0:                        
                        precomp.makeSynapse(postcomp, 
                                            name='ampa_from_%s' % (pretype_vertex['label']), 
                                            classname=synchan_classname, 
                                            Ek=0.0, 
                                            Gbar=g_ampa, 
                                            tau1=syn_edge['tauampa'], 
                                            tau2=syn_edge['tauampa'], 
                                            Pr=p_release, 
                                            delay=delay)
                    g_nmda = self.g_nmda_mat[pre_index, post_index]
                    if g_nmda != 0.0:
                        synchan = precomp.makeSynapse(postcomp, 
                                                      name='nmda_from_%s' % (pretype_vertex['label']), 
                                                      classname=nmdachan_classname, 
                                                      Ek=0.0, 
                                                      tau1=syn_edge['taunmda'], 
                                                      tau2=5e-3, 
                                                      Pr=p_release, 
                                                      delay=delay)
                        synchan.MgConc = TraubFullNetData.MgConc
                    g_gaba = self.g_gaba_mat[pre_index, post_index]
                    if g_gaba != 0.0:
                        g_gaba_slow = 0.0
                        if syn_edge['taugabaslow'] > 0.0:
                            if pretype_vertex['label'] == 'nRT':
                                if posttype_vertex['label'] == 'nRT':                                
                                    g_gaba_slow = g_gaba * (1 - self.frac_nRT_nRT_ggaba_fast)
                                    g_gaba = g_gaba * self.frac_nRT_nRT_ggaba_fast
                                elif posttype_vertex['label'] == 'TCR':
                                    g_gaba_slow = g_gaba *  (1 - self.frac_nRT_TCR_ggaba_fast)
                                    g_gaba = g_gaba * self.frac_nRT_TCR_ggaba_fast
                                precomp.makeSynapse(postcomp, 
                                                    name='gaba_slow_from_%s' % (pretype_vertex['label']), 
                                                    classname=synchan_classname, 
                                                    Ek=syn_edge['ekgaba'], 
                                                    Gbar=g_gaba_slow, 
                                                    tau1=syn_edge['taugabaslow'], 
                                                    tau2=0.0, 
                                                    Pr=p_release, 
                                                    delay=delay)                               

                        precomp.makeSynapse(postcomp, 
                                            name='gaba_from_%s' % (pretype_vertex['label']), 
                                            classname=synchan_classname, 
                                            Ek=syn_edge['ekgaba'], 
                                            Gbar=g_gaba,
                                            tau1=syn_edge['taugaba'], 
                                            tau2=0.0, 
                                            Pr=p_release, 
                                            delay=delay)
        endtime = datetime.now()
        delta = endtime - starttime
        config.BENCHMARK_LOGGER.info('Finished network creation in: %g s' % (delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))

    def setup_ectopic_input(self):
        config.LOGGER.debug('Setting up ectopic input')
        for celltype in self.celltype_graph.vs:
            for ii in self.populations[celltype['label']]:
                cell = self.index_cell_map[ii]
                randspike = moose.RandomSpike('ectopic_%s' % (cell.name), self.ectopic_container)
                randspike.rate = 1/celltype['ectopicinterval']
                randspike.minAmp = 0.4e-9
                randspike.maxAmp = 0.4e-9
                randspike.reset = 1
                randspike.resetValue = 0.0
                success = randspike.connect('outputSrc', cell.comp[cell.presyn], 'injectMsg')
                config.LOGGER.debug('Connected %s to %s: %s' % (randspike.path, cell.comp[cell.presyn].path, str(success)))

    def setup_lfp_recording(self, name, depth, data_container):
        """Setup electrodes for recording LFP."""
        electrode = moose.Efield(name, self.electrode_container)
        self.electrodes.append(electrode)
        electrode.x = 0.0
        electrode.y = 0.0
        electrode.z = depth
        for celltype in self.celltype_graph.vs:
            cellclass = eval(celltype['label'])
            if not hasattr(cellclass, 'depth') or cellclass.depth is  None:
                continue
            # first collect the indices of all the compartments that affect lfp
            comp_indices = []
            for level_no in cellclass.depth.keys():
                print
                print level_no, '-->',
                for comp_no in cellclass.level[level_no]:
                    print comp_no, ',',
                    comp_indices.append(comp_no)
            # then connect all such compartments in all cells of this type to the electrode object.
            for ii in self.populations[celltype['label']]:
                cell = self.index_cell_map[ii]
                # print 'cell:', cell.path, 'comps:', comp_indices
                for jj in comp_indices:
                    comp = cell.comp[jj]
                    result = comp.connect('ImSrc', electrode, 'currentDest')
                    config.LOGGER.debug('Connected %s %s [pos=(%g, %g, %g)] to %s %s: %s' % (comp.className, comp.path, comp.x, comp.y, comp.z, electrode.className, electrode.path, str(result)))

        lfp_container = moose.Neutral('lfp', data_container)
        lfp_table = moose.Table(name, lfp_container)
        lfp_table.stepMode = 3
        electrode.connect('potential', lfp_table, 'inputRequest')
        config.LOGGER.debug('Created electrode: %s at depth %g m' % (name, depth))
                               

    def get_maxdegree_cell_indices(self, celltype=None, size=None):
        """Get the cells with maximum connectivity - disregarding the strength of the synapse.

        returns the {size} cells of type {celltype} sorted by degree.
        
        """
        cell_dict = defaultdict(int)
        if celltype is not None:
            index = 0
            for vertex in self.celltype_graph.vs:
                if vertex['label'] == celltype:
                    for ii in range(index, index + vertex['count']):
                        for jj in range(self.g_gaba_mat.shape[1]):
                            if self.g_gaba_mat[ii, jj] != 0.0:
                                cell_dict[ii] += 1
                            if self.g_ampa_mat[ii, jj] != 0.0:
                                cell_dict[ii] += 1
                            if self.g_nmda_mat[ii, jj] != 0.0:
                                cell_dict[ii] += 1
                    break
                else:
                    index += vertex['count']
        else:
            for ii in range(self.g_gaba_mat.shape[0]):
                for jj in range(self.g_gaba_mat.shape[1]):
                    if self.g_gaba_mat[ii, jj] != 0.0:
                        cell_dict[ii] += 1
                    if self.g_ampa_mat[ii, jj] != 0.0:
                        cell_dict[ii] += 1
                    if self.g_nmda_mat[ii, jj] != 0.0:
                        cell_dict[ii] += 1
        cells = sorted(cell_dict, key=lambda x: cell_dict[x], reverse=True)
        if size is None:
            return cells
        else:
            return cells[:size]
            
    def setup_spike_recording(self, data_container):
        """Create tables to record spike times for each cell.

        The tables are created under {data_container}/spikes"""
        spike_container = moose.Neutral('spikes', data_container)
        for cell in self.cell_index_map.keys():
            tab = cell.soma.insertRecorder(cell.name, 'Vm', spike_container)
            tab.stepMode = moose.TAB_SPIKE
            tab.stepSize = 0.0
        for ch in self.ectopic_container.children():
            spike = moose.Neutral(ch)
            tab = moose.Table(spike.name, spike_container)
            tab.stepMode = moose.TAB_SPIKE
            tab.stepSize = 0.2e-9
            spike.connect('state', tab, 'inputRequest')
        
    def setup_Vm_recording(self, data_container, celltype, numcells=10, random=True):
        """Create tables to record Vm and [Ca2+] from a fixed number of cells of each type.

        The tables recording Vm are created under {data_container}/Vm and
        those for recording [Ca2+] are created under {data_container}/Ca.

        celltype -- the class of cells to record from.
        
        numcells -- number of cells to record from. If 'all', select all cells.

        random -- if true, randomly select {numcells} cells. if false, select only the maxdegree cells.    
        
        """
        vm_container = moose.Neutral('Vm', data_container)
        ca_container = moose.Neutral('Ca', data_container)
        if celltype == 'all':
            vs = self.celltype_graph.vs
        else:
            vs = self.celltype_graph.vs.select(label_eq=celltype)
        for vertex in vs:
            if not random:
                cell_list = self.get_maxdegree_cell_indices(celltype=vertex['label'], size=numcellspertype)
            else:
                pop = numpy.array(self.populations[celltype])
                high = len(pop)
                indices = numpy.random.randint(low=0, high=high, size=numcells)
                cell_list = pop[indices]
            for cellindex in cell_list:
                cell = self.index_cell_map[cellindex]
                cell.soma.insertRecorder(cell.name, 'Vm', vm_container)
                cell.soma.insertCaRecorder(cell.name, ca_container)
        
    def setup_stimulus(self, bg_celltype='any', 
                       probe_celltype='any',
                       stim_onset=1.0, 
                       bg_interval=0.5, 
                       pulse_width=60e-6, 
                       isi=10e-3, 
                       level=5e-12, 
                       bg_count=100, 
                       probe_count=10, 
                       stim_container='/stim'):
        """Setup the stimulus protocol.

        The protocol is as follows:

     ->|  |<- stim_onset
           __________________________________ gate
        __|
        ->|    |<- bg_interval

        _______||____||____||____||____||____  background
        
        _____________||__________||__________  probe


        Let the system stabilize for stim_onset seconds.

        Then turn the stim_gate on: which gates the triggers.

        The bg_trigger will trigger the background pulse
        generator. probe_trigger will trigger the probe pulse
        generator.

        stim_container -- container object for stimulating-electrodes

        celltype -- type of cells we are looking at

        stim_onset -- when we consider the system stabilized and start
        stimulus

        bg_interval -- interval between two stimulus sessions

        pulse_width -- width of background pulses

        isi -- if paired pulse, then the interval between the two
        pulses (beginning of second - beginning of first), 0 for
        single pulse.

        level -- current injection value.

        bg_count -- number of cells stimulated by background pulse.

        probe_count -- number of cells to be probed.

        """
        if  isinstance(stim_container, str):
            self.stim_container = moose.Neutral(stim_container)
        elif isinstance(stim_container, moose.Neutral):
            self.stim_container = stim_container
        else:
            raise Exception('Stimulus container must be a string or a Neutral object: got %s', stim_container.__class__.__name__)
        self.stim_gate = moose.PulseGen('stim_gate', self.stim_container)
        self.stim_gate.trigMode = moose.FREE_RUN
        self.stim_gate.firstDelay = stim_onset
        self.stim_gate.firstWidth = 1e9 # Keep it on forever
        self.stim_gate.firstLevel = 1.0                
        self.stim_bg = moose.PulseGen('stim_bg', self.stim_container)
        self.stim_bg.firstLevel = level
        self.stim_bg.secondLevel = level
        self.stim_bg.firstDelay = bg_interval
        self.stim_bg.firstWidth = pulse_width
        self.stim_bg.secondDelay = isi
        self.stim_bg.secondWidth = pulse_width        
        self.stim_bg.trigMode = moose.EXT_GATE

        self.stim_probe = moose.PulseGen('stim_probe', self.stim_container)
        self.stim_probe.firstLevel = level
        self.stim_probe.secondLevel = level
        self.stim_probe.firstDelay = 2 * bg_interval + pulse_width + isi
        self.stim_probe.secondDelay = isi
        self.stim_probe.firstWidth = pulse_width
        self.stim_probe.secondWidth = pulse_width            
        self.stim_probe.trigMode = moose.EXT_GATE
        self.stim_gate.connect('outputSrc', self.stim_bg, 'input')
        self.stim_gate.connect('outputSrc', self.stim_probe, 'input')
        bg_cell_indices = []
        probe_cell_indices = []
        if bg_celltype == 'any':
            bg_cell_indices = numpy.random.randint(low=0, high=len(self.index_cell_map.keys()), size=bg_count)
        else:
            celltype_vertex_set = self.celltype_graph.vs.select(label_eq=bg_celltype)
            for vertex in celltype_vertex_set:
                startindex = int(vertex['startindex'])
                count = int(vertex['count'])
                np.concatenate((bg_cell_indices, numpy.random.randint(low=startindex, high=startindex+count, size=bg_count)))

        probe_cell_indices = []
        if probe_celltype == 'any':
            probe_cell_indices = numpy.random.randint(low=0, high=len(self.index_cell_map.keys()), size=probe_count)
        else:
            celltype_vertex_set = self.celltype_graph.vs.select(label_eq=probe_celltype)
            for vertex in celltype_vertex_set:
                startindex = int(vertex['startindex'])
                count = int(vertex['count'])
                np.concatenate((bg_cell_indices, numpy.random.randint(low=startindex, high=startindex+count, size=probe_count)))
        protocol_file = open('protocol_'+ config.filename_suffix)
        protocol_file.write('onset: %g\nbg_interval: %g\nisi: %g\nwidth: %g\nlevel: %g\n' % (stim_onset, bg_interval, isi, pulse_width, level))
        protocol_file.write('background_cells:\n')
        for index in bg_cell_indices:            
            cell = self.index_cell_map[index]
            self.stim_bg.connect('outputSrc', cell, 'injectMsg')
            protocol_file.write('%s\n' % (cell.path))
        protocol_file.write('probe_cells:\n')
        for index in probe_cell_indices:
            cell = self.index_cell_map[index]
            self.stim_probe.connect('outputSrc', cell, 'injectMsg')
            protocol_file.write('%s\n' % (cell.path))
        protocol_file.close()

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
                pretype_vertex_seq = self.celltype_graph.vs
            else:
                pretype_vertex_seq = self.celltype_graph.vs.select(label_eq=celltype_pair[0])
            if celltype_pair[1] == '*':
                posttype_vertex_seq = self.celltype_graph.vs
            else:
                posttype_vertex_seq = self.celltype_graph.vs.select(label_eq=celltype_pair[1])
            for pretype_vertex in pretype_vertex_seq:
                for posttype_vertex in posttype_vertex_seq:
                    try:
                        edge_id = self.celltype_graph.get_eid(pretype_vertex.index, posttype_vertex.index)
                        edge = self.celltype_graph.es[edge_id]
                        if conductance_name in edge.attribute_names():
                            if conductance_name == 'ggaba' and pretype_vertex['label'] == 'nRT' and posttype_vertex['label'] == 'TCR':
                                self.nRT_TCR_ggaba_low *= scale_factor
                                self.nRT_TCR_ggaba_high *= scale_factor
                                edge[conductance_name] = 'uniform %f %f' % (self.nRT_TCR_ggaba_low, self.nRT_TCR_ggaba_high)
                            else:
                                edge[conductance_name] *= scale_factor

                    except Exception, e:
                        pass

    def set_populations(self, filename):
        """Read cellcounts from specified file and updates the
        celltype-graph accordingly. The file should have space
        separated values like:

        name count

        """
        if filename is None:
            return
        with open(filename) as popcount_file:
            config.LOGGER.info('Reading cell counts from: %s' % (filename))
            for line in popcount_file.readlines():
                if line.strip().startswith('#'):
                    continue
                tokens = line.split()
                if not tokens:
                    continue
                cellname, count = tokens
                vertices = self.celltype_graph.vs.select(label_eq=cellname)
                for vertex in vertices:
                    vertex['count'] = int(count)
                    config.LOGGER.info('%s population size: %d' % (cellname, vertex['count']))

    def tweak_Ek(self, channel_class, value):
        """Adds value to channel's reversal potential. According
        Nernst equation, a multiplicative change in ionic
        concentration will cause an additive change in the reversal
        potential for that ion."""
        config.LOGGER.info('%s Ek += %g' % (channel_class.__name__, value))
        self.tweaks_doc.append('%s.Ek += %g' % (channel_class.__name__, value))
        for cell in self.cell_index_map.keys():
	    for ii in range(1, cell.num_comp):
                comp = cell.comp[ii]
                for chan in comp.channels:
                    if isinstance(chan, channel_class):
                        chan.Ek += value

    def set_unknown_prelease(self, value):
        self.tweaks_doc.append('prelease = 1.0 <- %g' % (value))
        for syn in self.celltype_graph.es:
            if syn['prelease'] == 1.0:
                syn['prelease'] = value

    def tune_conductances(self, filename):
        """Tune the conductances based on entries in file filename.
        file should have space separated entries like this:

        conductance-name sourcetype desttype scalefactor
        """
        if filename is None:
            return
        with open(filename) as synfile:
            for line in synfile.readlines():
                if line.strip().startswith('#'):
                    continue
                tokens = line.split()
                if not tokens:
                    continue
                [g_name, source, dest, scale_factor] = tokens                
                scale_factor = float(scale_factor)
                source_vertices = self.celltype_graph.select(label_eq=source)
                dest_vertices = self.celltype_graph.select(label_eq=dest)
                for src_v in source_vertices:
                    for dest_v in dest_vertices:
                        synid = self.celltype_graph.get_eid(src_v, dest_v)
                        if synid:
                            syn = self.celltype_graph.es[synid]
                            if g_name in syn.attribute_names():
                                if g_name == 'ggaba' and source == 'nRT' and dest == 'TCR':
                                    self.nRT_TCR_ggaba_low *= scale_factor
                                    self.nRT_TCR_ggaba_high *= scale_factor
                                    syn[g_name] = 'uniform %f %f' % (self.nRT_TCR_ggaba_low, self.nRT_TCR_ggaba_high)

                                else:
                                    syn[g_name] *= scale_factor
                                self.tweaks_doc.append('%s[%s->%s] *= %g' % (source, dest, scale_factor))

                                    
    def set_conductances(self, filename):
        """Set the conductances based on entries in file filename.
        file should have space separated entries like this:

        conductance-name sourcetype desttype value
        """
        if filename is None:
            return
        with open(filename) as synfile:
            for line in synfile.readlines():
                if line.strip().startswith('#'):
                    continue
                tokens = line.split()
                if not tokens:
                    continue
                [g_name, source, dest, value] = tokens
                value = float(value)
                source_vertices = self.celltype_graph.select(label_eq=source)
                dest_vertices = self.celltype_graph.select(label_eq=dest)
                for src_v in source_vertices:
                    for dest_v in dest_vertices:
                        synid = self.celltype_graph.get_eid(src_v, dest_v)
                        if synid:
                            syn = self.celltype_graph.es[synid]
                            if g_name in syn.attribute_names():
                                if g_name == 'ggaba' and source == 'nRT' and dest == 'TCR':
                                    self.nRT_TCR_ggaba_low = value
                                    self.nRT_TCR_ggaba_high = value
                                    syn[g_name] = 'uniform %f %f' % (self.nRT_TCR_ggaba_low, self.nRT_TCR_ggaba_high)
                                else:
                                    syn[g_name] = value


    def save_network_model(self,  filename):
        """Save the network structure in an hdf5 file"""
        config.LOGGER.debug('Start saving the network model')
        starttime =  datetime.now()
        compression_filter =  tables.Filters(complevel=9, complib='zlib', fletcher32=True)
        h5file =  tables.openFile(filename,  mode = 'w',  title = 'Traub Network: timestamp: %s' % (config.timestamp.strftime('%Y-%M-%D %H:%M:%S')),  filters = compression_filter)
        h5file.root._v_attrs.rngseed = config.rngseed
        h5file.root._v_attrs.notes = '\n'.join(self.tweaks_doc)
        # Save the celltype information (vertices of the celltype graph)
        network_struct =  h5file.createGroup(h5file.root, 'network', 'Network structure')
        celltype_table =  h5file.createTable(network_struct, 'celltype', CellType,  'Information on each celltype population')
        celltype =  celltype_table.row
        for vertex in self.celltype_graph.vs:
            celltype['name'] =  vertex['label']
            celltype['index'] =  vertex.index
            celltype['count'] =  vertex['count']
            celltype.append()
        synapse_table = h5file.createTable(network_struct,  'synapsetype',  SynEdge, 'Synapse information between celltype pairs')
        synedge =  synapse_table.row
        for edge in self.celltype_graph.es:
            synedge['source'] = edge.source
            synedge['target'] =  edge.target
            synedge['weight'] = edge['weight']
            synedge['gampa'] = edge['gampa']
            synedge['gnmda'] = edge['gnmda']
            synedge['tauampa'] = edge['tauampa']
            synedge['taunmda'] = edge['taunmda']
            synedge['tau2nmda'] =  5e-3
            synedge['taugaba'] = edge['taugaba']
            synedge['prelease'] = edge['prelease']
            ii =  0
            pscomps = numpy.zeros(90, dtype=numpy.uint8)
            for pscomp in eval(edge['pscomps']): 
                pscomps[ii] = int(pscomp)
                ii +=  1
            synedge['pscomps'] = pscomps
            synedge['ekgaba'] = edge['ekgaba']
            # print edge.source, edge.target, edge['ggaba'], type(edge['ggaba'])

            it =  None
            try:
                it =  iter(edge['ggaba'])
            except TypeError:
                synedge['ggaba'] = numpy.array([edge['ggaba'], edge['ggaba']])

            assert ((it is None) or (self.celltype_graph.vs[edge.source]['label'] == 'nRT'))
            if self.celltype_graph.vs[edge.source]['label'] == 'nRT':
                if self.celltype_graph.vs[edge.target]['label'] == 'TCR':
                    synedge['ggaba'] =  numpy.array([self.nRT_TCR_ggaba_low, self.nRT_TCR_ggaba_high])
                synedge['taugabaslow'] = edge['taugabaslow']
            synedge.append()
        cellnet_group = h5file.createGroup(network_struct, 'cellnetwork', 'Cell-to-cell network structure')
        if self.g_ampa_mat.nnz > 0:
            gampa_array =  h5file.createCArray(cellnet_group, 'gampa', tables.FloatAtom(),  shape=(self.g_ampa_mat.nnz, 3))
            ii =  0
            for (index,  value) in self.g_ampa_mat.items():
                gampa_array[ii,0] = index[0]
                gampa_array[ii,1] =  index[1]
                gampa_array[ii,2] = value
                ii +=  1
        if self.g_nmda_mat.nnz > 0:
            gnmda_array =  h5file.createCArray(cellnet_group, 'gnmda', tables.FloatAtom(),  shape=(self.g_nmda_mat.nnz, 3))
            ii =  0
            for (index,  value) in self.g_nmda_mat.items():
                gnmda_array[ii,0] = index[0]
                gnmda_array[ii,1] =  index[1]
                gnmda_array[ii,2] = value
                ii +=  1
        if self.g_gaba_mat.nnz > 0:
            ggaba_array =  h5file.createCArray(cellnet_group, 'ggaba', tables.FloatAtom(),  shape=(self.g_gaba_mat.nnz, 3))
            ii =  0
            for (index,  value) in self.g_gaba_mat.items():
                ggaba_array[ii,0] = index[0]
                ggaba_array[ii,1] =  index[1]
                ggaba_array[ii,2] = value
                ii +=  1
        if self.ps_comp_mat.nnz > 0:
            pscomp_array =  h5file.createCArray(cellnet_group, 'pscomp',  tables.Int32Atom(),  shape=(self.ps_comp_mat.nnz, 3))
            ii =  0
            for (index,  value) in self.ps_comp_mat.items():
                pscomp_array[ii,0] = index[0]
                pscomp_array[ii,1] =  index[1]
                pscomp_array[ii,2] = value
                ii +=  1
        h5file.close()
        endtime =  datetime.now()
        delta =  endtime -  starttime
        config.BENCHMARK_LOGGER.info('Saved network model in:% g s' %  (delta.days *  86400 +  delta.seconds +  1e-6 * delta.microseconds))

    def verify_saved_model(self, filename):
        starttime = datetime.now()
        h5file =  tables.openFile(filename)
        celltypes =  h5file.getNode('/network', name='celltype')
        for row in celltypes.iterrows():
            index =  row['index']
            assert self.celltype_graph.vs[index]['label'] == row['name']
            assert self.celltype_graph.vs[index]['count'] == row['count']
        synedges =  h5file.getNode('/network', 'synapsetype')
        for row in synedges.iterrows():
            source = row['source']
            target =  row['target']
            edge = self.celltype_graph.es[self.celltype_graph.get_eid(source, target)]
            assert row['ekgaba'] == edge['ekgaba']
            assert row['weight'] == edge['weight']
            assert row['gampa'] == edge['gampa']
            assert row['gnmda'] == edge['gnmda']
            assert row['tauampa'] == edge['tauampa']
            assert row['taunmda'] == edge['taunmda']
            assert row['tau2nmda'] ==  5e-3
            assert row['taugaba'] == edge['taugaba']
            ii =  0
            for pscomp in eval(edge['pscomps']): 
                try:
                    assert row['pscomps'][ii] == pscomp
                except AssertionError:
                    config.LOGGER.debug('pscomp not same: %s <> %s in synapse from %s to %s' % (row['pscomps'][ii], pscomp, self.celltype_graph.vs[source]['label'], self.celltype_graph.vs[target]))
                ii +=  1
            assert row['ekgaba'] == edge['ekgaba']

            it =  None
            try:
                it =  iter(edge['ggaba'])
            except TypeError:
                assert row['ggaba'][0] == edge['ggaba']
                assert row['ggaba'][0] == edge['ggaba']

            assert ((it is None) or (self.celltype_graph.vs[edge.source]['label'] == 'nRT'))
            if self.celltype_graph.vs[edge.source]['label'] == 'nRT':
                if self.celltype_graph.vs[edge.target]['label'] == 'TCR':
                    assert row['ggaba'][0] == self.nRT_TCR_ggaba_low
                    assert row['ggaba'][1] == self.nRT_TCR_ggaba_high
                assert row['taugabaslow'] == edge['taugabaslow']
        h5file.close()
        endtime = datetime.now()
        delta = endtime - starttime
        config.BENCHMARK_LOGGER.info('Finished verification of saved model in hdf5 in: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        config.LOGGER.info('Verified model in: %s :: SUCCESS' %(filename))

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


def test_reading_network(filename):
    tn =  TraubNet()
    tn._generate_celltype_graph()
    tn._generate_cell_graph()
    h5file =  tables.openFile(filename)
    celltypes =  h5file.getNode('/network', name='celltype')
    for row in celltypes.iterrows():
        index =  row['index']
        assert tn.celltype_graph.vs[index]['label'] == row['name']
        assert tn.celltype_graph.vs[index]['count'] == row['count']
    synedges =  h5file.getNode('/network', 'synapsetype')
    for row in synedges.iterrows():
        source = row['source']
        target =  row['target']
        egde = tn.celltype_graph.es[tn.get_eid(source, target)]
        assert edge['ekgaba'] == edge['ekgaba']
        assert row['weight'] == edge['weight']
        assert row['gampa'] == edge['gampa']
        assert row['gnmda'] == edge['gnmda']
        assert row['tauampa'] == edge['tauampa']
        assert row['taunmda'] == edge['taunmda']
        assert row['tau2nmda'] ==  5e-3
        assert row['taugaba'] == edge['taugaba']
        ii =  0
        for pscomp in eval(edge['pscomps']): 
            assert row['pscomps'][ii] == pscomp
            ii +=  1
        assert row['ekgaba'] == edge['ekgaba']

        it =  None
        try:
            it =  iter(edge['ggaba'])
        except TypeError:
            assert row['ggaba'][0] == edge['ggaba']
            assert row['ggaba'][0] == edge['ggaba']

        assert ((it is None) or (self.celltype_graph.vs[edge.source]['label'] == 'nRT'))
        if self.celltype_graph.vs[edge.source]['label'] == 'nRT':
            if self.celltype_graph.vs[edge.target]['label'] == 'TCR':
                assert row['ggaba'][0] ==  tn.nRT_TCR_ggaba_low
                assert row['ggaba'][1] ==  tn.nRT_TCR_ggaba_high
            assert row['taugabaslow'] == edge['taugabaslow']
        
    h5file.close()

if __name__ == '__main__':
    net = TraubNet()
    net._generate_celltype_graph()
    net._generate_cell_graph()
    net.create_network()
    net.save_network_model(config.MODEL_FILENAME)
    net.verify_saved_model(config.MODEL_FILENAME)
    
# 
# trbnet.py ends here
