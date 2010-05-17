# network.py --- 
# 
# Filename: network.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
# Last-Updated: Mon May 17 15:36:59 2010 (+0530)
#           By: subha
#     Update #: 515
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This is for setting up the network. The connectivity is
# stored as a map of maps.
# 
# 
# 
# 

# Change log:
#
# 2010-01-13 initial version to load csv
#
# 2010-01-14 added cell count for the original model
# 
# 2010-03-29 12:15:20 (+0530) -- shifted the original file to network_spec.py
#                                the current file will be used for creating 
#                                the complete network - more as a driver program.

# Code:

from numpy import *
from datetime import datetime

import moose
import config

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

from simulation import Simulation
from population import Population

CELL_COUNT = {
    'SupPyrRS': 1000,
    'SupPyrFRB': 50,
    'SupBasket': 90,
    'SupAxoaxonic': 90,
    'SupLTS': 90,
    'SpinyStellate': 240,
    'TuftedIB': 800,
    'TuftedRS': 200,
    'NontuftedRS': 500,
    'DeepBasket': 100,
    'DeepAxoaxonic': 100,
    'DeepLTS': 100,
    'TCR': 100,
    'nRT': 100
}

class Network:
	"""A network of cell-populations of different kinds"""
	def __init__(simulation, pop_counts):
		"""init method of the network class.

		simulation -- simulation object giving handle to the containers.

		pops -- a dictionary containing each cell class and its population size.

		This is a rather specific way of designing it. Each
		population will have the same name as the cell class.		

		"""
		self.population_count = pop_counts
		self.sim = simulation
		self.population_list = []
		start = datetime.now()
 
		for cell_type, count in pop_counts.items():
			cell_class = eval(cell_type)
			self.population_list.append(Population(sim.model.path + '/' + cell_type, cell_class, count))
            for pre_population in self.population_list:
				for post_population in self.population_list:
					pre_population.connect(post_population)
 
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('time to create all the population: %g' % (delta.seconds + 1e-6 * delta.microseconds))


    def setup_random_recording(self, n=1):
        """Setup a bunch of tables to record from random cells in each population.
        
        Select some cells from each population and record the Vm.
        
        n -- number/proportion of cells to be recorded from. If n is
        an integer then it is the absolute number of cells from which
        the Vm will be recorded. If it is a float <= 1 in abosulte
        value, it represents the proportion of cells from which the
        recording should be done.

        """
        
        if isinstance(n, float) and abs(n) <= 1.0:
            for pop in self.population_list:
                presyn = pop.cell_class.presyn
                count = n * len(pop.cell_list)
                cell_no_list = random.randint(count)
                for cell_no in cell_no_list:
                    cell = pop.cell_list[cell_no]
                    comp = cell.comp[presyn]
                    recorder_name = '%s/%s__%d' % (self.sim.data.path, cell.name, presyn)
                    recorder = moose.Table(recorder_name)
                    recorder.stepMode = 3
                    recorder.connect('inputRequest', comp, 'Vm')
        elif isinstance(n, int):
            for pop in self.population_list:
                presyn = pop.cell_class.presyn
                count = n
                cell_no_list = random.randint(0, high=count, size=count)
                for cell_no in cell_no_list:
                    cell = pop.cell_list[cell_no]
                    comp = cell.comp[presyn]
                    recorder_name = '%s/%s__%d' % (self.sim.data.path, cell.name, presyn)
                    recorder = moose.Table(recorder_name)
                    recorder.stepMode = 3
                    recorder.connect('inputRequest', comp, 'Vm')

            

def test_full_model(simtime, simdt=1e-4, plotdt=1e-3):
    """Setup and run the full Traub model"""
    sim = Simulation('traub')
    net = []
    scale = 10
    
    network = Network(sim, CELL_COUNT)
    network.setup_random_recording(0.05) # record from 5% of cells
    sim.schedule(simdt=simdt, plotdt=plotdt, gldt=1e10)
    sim.run(time=simtime)
    sim.dump_data('data', True)

def test_all_cell_type():
    """test-load all different cell type. this is for debugging - as test_full_model is crashing silently after reading nRT cell"""
    cells = []
    for cell_type in CELL_COUNT.keys():
        print '####', cell_type
        cell_class = eval(cell_type)
	cells.append(cell_class(cell_class.prototype, cell_type))
	config.LOGGER.debug('Created cell %s' % (cell_type))


if __name__ == '__main__':
    # test_all_cell_type()
    test_full_model(1000e-3)
    config.LOGGER.info('Finished simulation')


    
# 
# network.py ends here
