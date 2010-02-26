# population.py --- 
# 
# Filename: population.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Feb 18 22:00:46 2010 (+0530)
# Version: 
# Last-Updated: Sat Feb 27 00:05:36 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 477
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# 
# 
# 

# Change log:
# 
# 
# 
# 

# Code:

from collections import defaultdict
from datetime import datetime
import csv 
import numpy
import moose
import allowedcomp
import synapse
import config
from simulation import Simulation
class Population(moose.Neutral):
    """Homogeneous cell population - handles setting up connections
    between populations based on connection matrix.

    cell_class -- class name of cells contained in this population
    cell_list -- cells contained in this population
    conn_map -- a dictionary of actual connection map for each post-synaptic population.
                the key is the target Population object, value is an n X m x 2 array where n
                is the number of cells in the postsynaptic population, m is the number of 
                presynaptic cells from this population per postsynaptic cell in the target 
                population.

                conn_map[target][i][j][0] is
                for the i-th cell in 'target' population, the index of the j-th presynaptic 
                cell in this population.

                conn_map[target][i][j][1] is 
                for the i-th cell in 'target' population, the index of the postsynaptic 
                compartment for j-th presynaptic cell in this population.

                This is rather convoluted - but there is no clear workaround because of the
                way the connectivity is specified (number of pre-cells per post-cell).

    """
    

    CELL_CONNECTION_MAP = None
    ALLOWED_COMP_MAP = None

    def __init__(self, path, cell_class, cell_count, prefix=None):
	"""Initialze the population by creating the cells.

	path -- path of this element as a container
	cell_class -- classobject for the cell type to be contained in
	this element.
	cell_count -- number of cells in this population
	prefix -- the n-th cell in this population gets the name {prefix}_{n}

	"""
	moose.Neutral.__init__(self, path)
        self.get_allowed_comp_map()
        self.get_connection_map()
	self.cell_list = []
	self.cell_type = cell_class.__name__
        self.cell_class = cell_class
	if prefix is None:
	    prefix = self.cell_type
	for number in range(cell_count):
	    cell_name = '%s_%d' % (prefix, number)
	    cell_instance = cell_class(cell_class.prototype, self.path + '/' + cell_name)
	    self.cell_list.append(cell_instance)
        self.conn_map = {}
        self.glView = None

    def connect(self, target):
        """Connect cells from this population to cells on the target
        population.        
        
        """
        config.LOGGER.debug(__name__ + ' starting')
        start = datetime.now()
        connection_map = self.get_connection_map()
        num_pre_per_post = connection_map[self.cell_type][target.cell_type]
        self.conn_map[target] = numpy.zeros((len(target.cell_list), num_pre_per_post, 2), dtype=int)
        # This gets a 2_D matrix whose row[i] is the array of indices
        # of the pre-synaptic cells for i-th post-synaptic cell.
        precell_indices = numpy.random.randint(0, 
                                               high=len(self.cell_list), 
                                               size=(len(target.cell_list), num_pre_per_post))

        allowed_comp_map = self.get_allowed_comp_map()
        allowed_comp_list = numpy.array(allowed_comp_map[self.cell_type][target.cell_type], dtype=int)
        target_comp_indices = numpy.random.randint(0, 
                                                   high=len(allowed_comp_list), 
                                                   size=(len(target.cell_list), num_pre_per_post))
        tau_GABA_fast = self.get_tauGABA(target)
        tau_GABA_slow = self.get_tauGABA(target, fast=False)
        tau_NMDA = self.get_tauNMDA(target)
        tau_AMPA = self.get_tauAMPA(target)
        gbar_GABA = self.get_GbarGABA(target)
        gbar_AMPA = self.get_GbarAMPA(target)
        gbar_NMDA = self.get_GbarNMDA(target)

        for ii in range(len(target.cell_list)):
            target_comp_list = allowed_comp_list[target_comp_indices[ii]] # list containing the target compartment for each presynaptic cell
            for jj in range(num_pre_per_post):
                precell_index = precell_indices[ii][jj]
                precell = self.cell_list[precell_index]
                postcell = target.cell_list[ii]
                post_syn_comp_index = target_comp_list[jj]
                post_syn_comp = postcell.comp[post_syn_comp_index]
                pre_syn_comp = precell.comp[precell.presyn]
                config.LOGGER.debug('connecting: \t%s \tto \t%s' % (pre_syn_comp.path, post_syn_comp.path))

                if tau_GABA_fast is not None:
                    pre_syn_comp.makeSynapse(post_syn_comp, name='GABA', Ek=self.cell_class.EGABA, Gbar=gbar_GABA, tau1=tau_GABA_fast, tau2=tau_GABA_slow, absRefract=1.5e-3)
                    config.LOGGER.debug('%s\tto%s\tGABA' % (pre_syn_comp.path, post_syn_comp.path))
                if tau_AMPA is not None:
                    pre_syn_comp.makeSynapse(post_syn_comp, name='AMPA', Ek=0.0, tau1=tau_AMPA, tau2=tau_AMPA)
                    config.LOGGER.debug('%s\tto%s\tAMPA' % (pre_syn_comp.path, post_syn_comp.path))
                if tau_NMDA is not None:
                    # TODO: NMDA time course is defined as:
                    # c * g(V, [Mg2+]) * S(t) where c is scaling constant,
                    # g is a function of V and [Mg2+]o, 0 < g <= 1 
                    # and S(t) is time dependent component of the ligand gated conductance.
                    pre_syn_comp.makeSynapse(post_syn_comp, name='NMDA', Ek=0.0, tau1=tau_NMDA, tau2=tau_NMDA, absRefract=1.5e-3)
                    config.LOGGER.debug('%s\tto%s\tNMDA' % (pre_syn_comp.path, post_syn_comp.path))

                self.conn_map[target][ii][jj][0] = precell_index
                self.conn_map[target][ii][jj][1] = post_syn_comp_index
                
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('(%s[%d], %s[%d]) - time: %g' % (self.cell_type, len(self.cell_list), target.cell_type, len(target.cell_list), delta.seconds + 1e-6 * delta.microseconds))

    def get_connection_map(self, filename='connmatrix.txt'):
	"""Load the celltype-to-celltype connectivity map from file
	and return a nested dictionary dict where dict[X][Y] is the
	number of presynaptic cells of type X connecting to each
	postsynaptic cell of type Y.

	filename -- the path of a csv file containing the connectivity
	matrix. The first row in the file should be the column headers
	- which are the cell types. The connection matrix itself is a
	square matrix with entry[i][j] specifying the number of
	presynaptic cells per postsynaptic cell, where the presynaptic
	cells are of type header[i] and the postsynaptic cell is of
	type header[j]

	"""
        if Population.CELL_CONNECTION_MAP is not None:
            return Population.CELL_CONNECTION_MAP
	config.LOGGER.debug('load_connection_map - loading connection matrix from ' + filename)
	connection_map = defaultdict(dict)
	reader = csv.reader(file(filename))
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
		connection_map[pre][post] = value
		col += 1
        Population.CELL_CONNECTION_MAP = connection_map
	config.LOGGER.debug('load_connection_map - finished loading.')
	return connection_map
    
    def get_allowed_comp_map(self, filename=None):
	"""Load the tables for allowed compartment list for synapses
	between each types of cells. 
	
	Return a nested dictionary with entries like:
	entry["X"]["Y"] = [n1, n2, n3, n4, ...]

	where only compartments with index n1, n2, n3, .... are
	allowed as postsynaptic compartment for synapses from cells of
	type X to cells of type Y.

	filename -- the source of the map. A netcdf-4 file?

	"""
	# netCDF4 seems to be too complicated for simple things like a
        # map. So just using the dict objects
        if Population.ALLOWED_COMP_MAP is not None:
            return Population.ALLOWED_COMP_MAP
        if filename is None:
            Population.ALLOWED_COMP_MAP = allowedcomp.ALLOWED_COMP
        else:
            raise Error, '%s - Loading allowed_compartment map from file is not yet implemented' % (__name__)
        return Population.ALLOWED_COMP_MAP

    def get_tauGABA(self, post, fast=True):
        tau_gaba_map = None
        if self.cell_type != 'nRT':
            tau_gaba_map = synapse.TAU_GABA
        elif fast:
            tau_gaba_map = synapse.TAU_GABA_FAST
        else:
            tau_gaba_map = synapse.TAU_GABA_SLOW
        try:
            return tau_gaba_map[self.cell_type][post.cell_type]
        except KeyError:
            return None

    def get_tauAMPA(self, post):
        try:
            return synapse.TAU_AMPA[self.cell_type][post.cell_type]
        except KeyError:
            return None
    
    def get_tauNMDA(self, post):
        try:
            return synapse.TAU_NMDA[self.cell_type][post.cell_type]
        except KeyError:
            return None

    def get_GbarAMPA(self, post):
        try:
            return synapse.G_AMPA[self.cell_type][post.cell_type]
        except KeyError:
            return None
    
    def get_GbarNMDA(self, post):
        try:
            return synapse.G_NMDA[self.cell_type][post.cell_type]
        except KeyError:
            return None

    def get_GbarGABA(self, post):
        try:
            return synapse.G_GABA[self.cell_type][post.cell_type]
        except KeyError:
            return None
    
        

    def setup_visualization(self, glviewname, parent, host='localhost', port='9999'):
        self.glView = moose.GLview(glviewname, parent)
        self.glView.vizpath = self.path + '/#/comp_1'
        self.glView.port = port
        self.glView.host = host
        self.glView.value1 = 'Vm'
        self.glView.value1min = -0.1
        self.glView.value1max = 0.05
        self.glView.morph_val = 1
        self.glView.color_val = 1
        self.glView.sync = 'off'
        self.glView.grid = 'on'

from spinystellate import SpinyStellate
# from suppyrRS import SupPyrRS

import time
from glclient import GLClient
def start_test_client():
    testclient = GLClient(exe='/home/subha/src/moose/gl/src/glclient', mode='v', colormap='/home/subha/src/moose/gl/colormaps/rainbow2', save_directory='/tmp')
    return testclient

def test_main():
    client = start_test_client()
    # time.sleep(3)
    sim = Simulation('/sim')
    cellcount = 4
    start = datetime.now()
    pre = Population(sim.model.path + '/ss', SpinyStellate, cellcount)
    pre.setup_visualization('gl_' + pre.name, sim.data)
    end = datetime.now()
    delta = end - start
    config.BENCHMARK_LOGGER.info('time to create population of %d cells: %g' % (cellcount, delta.seconds + 1e-6 * delta.microseconds))
    
    post = pre
    pre.connect(post)
    precell_index = pre.conn_map[post][0][0][0]
    post_comp_index = pre.conn_map[post][0][0][1]
    print precell_index, post_comp_index
    precell = pre.cell_list[precell_index]
    precomp = precell.comp[precell.presyn]
    postcomp = post.cell_list[0].comp[post_comp_index]
    # precell.soma.insertPulseGen('inject', sim.model, firstLevel=0.0)
    preVmTable = precomp.insertRecorder('preVmTable', 'Vm', sim.data)
    postVmTable = postcomp.insertRecorder('postVmTable', 'Vm', sim.data)
    sim.schedule()
    sim.run()
    preVmTable.dumpFile('preVm.txt')
    postVmTable.dumpFile('postVmTable.txt')
    client.stop()

if __name__ == '__main__':
    test_main()
# 
# population.py ends here
