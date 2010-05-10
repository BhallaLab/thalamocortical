# network.py --- 
# 
# Filename: network.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed Jan 13 22:33:35 2010 (+0530)
# Version: 
# Last-Updated: Fri May  7 18:07:13 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 397
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

from datetime import datetime
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
              

def test_full_model(simtime, simdt=1e-4, plotdt=1e-3):
    """Setup and run the full Traub model"""
    sim = Simulation('traub')
    net = []
    scale = 10
    start = datetime.now()
 
    for cell_type, count in CELL_COUNT.items():
        cell_class = eval(cell_type)
        # net.append(Population(sim.model.path + '/' + cell_type, cell_class, count))
        net.append(Population(sim.model.path + '/' + cell_type, cell_class, 1))
    for pre_population in net:
        for post_population in net:
            pre_population.connect(post_population)
 
    end = datetime.now()
    delta = end - start
    config.LOGGER.info('time to create all the population: %g' % (delta.seconds + 1e-6 * delta.microseconds))
    sim.schedule(simdt=simdt, plotdt=plotdt, gldt=1e10)
    sim.run(time=simtime)

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
