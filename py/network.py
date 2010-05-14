# network.py --- 
# 
# Filename: network.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
<<<<<<< HEAD:py/network.py
# Last-Updated: Thu May 13 11:46:14 2010 (+0530)
#           By: subha
#     Update #: 482
=======
# Last-Updated: Mon May 10 11:28:58 2010 (+0530)
#           By: subha
#     Update #: 471
>>>>>>> 10e5f490fd90d554b63029bcb6e2eb2ddfcac36c:py/network.py
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

def setup_random_inject(simulation, cells_list, n=1):
	"""Set-up current injection."""
	pass

def setup_random_recording(simulation, pop_list, n=1):
        """Setup a bunch of tables to record from random cells in each population.

        Select some cells from each population and record the Vm.
        
        pop_list -- list of Population instances.
        
        n -- number/proportion of cells to be recorded from. If n is
        an integer then it is the absolute number of cells from which
        the Vm will be recorded. If it is a float <= 1 in abosulte
        value, it represents the proportion of cells from which the
        recording should be done.
        
        """
        
        if isinstance(n, float) and abs(n) <= 1.0:
                for pop in pop_list:
                        presyn = pop.cell_class.presyn
                        count = n * len(pop.cell_list)
                        cell_no_list = random.randint(count)
                        for cell_no in cell_no_list:
                                cell = pop.cell_list[cell_no]
                                comp = cell.comp[presyn]
                                recorder_name = '%s/%s__%d' % (simulation.data.path, cell.name, presyn)
                                recorder = moose.Table(recorder_name)
                                recorder.stepMode = 3
                                recorder.connect('inputRequest', comp, 'Vm')
        elif isinstance(n, int):
                for pop in pop_list:
                        presyn = pop.cell_class.presyn
                        count = n
                        cell_no_list = random.randint(0, high=count, size=count)
                        for cell_no in cell_no_list:
                                cell = pop.cell_list[cell_no]
                                comp = cell.comp[presyn]
                                recorder_name = '%s/%s__%d' % (simulation.data.path, cell.name, presyn)
                                recorder = moose.Table(recorder_name)
                                recorder.stepMode = 3
                                recorder.connect('inputRequest', comp, 'Vm')

            

def test_full_model(simtime, simdt=1e-4, plotdt=1e-3):
    """Setup and run the full Traub model"""
    sim = Simulation('traub')
    net = []
    scale = 10
    start = datetime.now()
 
    for cell_type, count in CELL_COUNT.items():
        cell_class = eval(cell_type)
        net.append(Population(sim.model.path + '/' + cell_type, cell_class, count))
    for pre_population in net:
        for post_population in net:
            pre_population.connect(post_population)
 
    end = datetime.now()
    delta = end - start
    config.LOGGER.info('time to create all the population: %g' % (delta.seconds + 1e-6 * delta.microseconds))
    setup_random_recording(sim, net)
    sim.schedule(simdt=simdt, plotdt=plotdt, gldt=1e10)
    sim.run(time=simtime)
    sim.dump_data('data')

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
    test_full_model(50e-3)
    config.LOGGER.info('Finished simulation')


    
# 
# network.py ends here
