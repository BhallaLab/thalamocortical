"""Simulate a specified cell with multiple current injection values
and plot the fI curve"""
import sys
from pylab import cm
from matplotlib import pyplot as plt
import numpy as np
import moose
import config
from simulation import Simulation
from spinystellate import SpinyStellate
from deepbasket import DeepBasket
from deepLTS import DeepLTS

cmap = cm.jet

def setup_model(testid, celltype, current, sim):
    """Setup a single cell model with specified current injection."""
    testid = '%s_%s' % (celltype, testid)
    cellclass = eval(celltype)
    cell = cellclass(cellclass.prototype, '%s/%s' % (sim.model.path, testid))
    config.LOGGER.debug('Created cell: %s' % (cell.path))
    vm_table = cell.comp[cell.presyn].insertRecorder('Vm_%s' % (testid), 'Vm', sim.data)
    # Create pulsegen with defaults.
    pulsegen = cell.soma.insertPulseGen('inject_%s' % (testid), sim.model, firstLevel=current, firstDelay=100e-3, firstWidth=1.0)
    im_table = cell.soma.insertRecorder('Im_%s' % (testid), 'Im', sim.data)
    return {'cell': cell,
            'Vm': vm_table,
            'Im': im_table,
            'inject': pulsegen
            }

def inject_test(celltype, current_list, simtime=2.0):
    test_data = {}
    sim = Simulation(celltype)
    for ii, current in enumerate(current_list):
        test_data[ii] = setup_model(ii, celltype, current, sim)
    sim.schedule()
    sim.run(simtime)
    sim.dump_data('data')
    timeseries = np.linspace(0, simtime, len(test_data[0]['Vm']))
    for ii, current in enumerate(current_list):
        color = cmap(ii*1.0/len(current_list), len(current_list))
        plt.plot(timeseries * 1e3, np.array(test_data[ii]['Vm']) * 1e3, color=color, label='%g nA' % (current * 1e9))
    plt.legend()
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print """Usage: %s celltype current0 [current1 [...]]

simulate current clamps of currentN nA on celltype."""
    celltype = sys.argv[1]
    currents = [float(ii)*1e-9 for ii in sys.argv[2:]]
    inject_test(celltype, currents)
