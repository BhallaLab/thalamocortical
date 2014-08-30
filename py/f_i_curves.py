"""obtain frequency vs current curve for single cells in traub model."""
import sys
import numpy as np
import h5py as h5

import config
import moose
from compartment import MyCompartment
from cell import TraubCell
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

DELAY = 500e-3
WIDTH = 500e-3
SIMTIME = 2.0
simdt = 0.25e-4
plotdt = 1e-4

model_container = moose.Neutral('/model')
data_container = moose.Neutral('/data')

def current_step_test(celltype, low, high, step, h5file):
    """Do current step test on cell type

    """
    cellclass = eval(celltype)
    rec_info = {}
    for ii in range(int((high - low) * 1.0 / step + 0.5)):
        test_id = '%s_%d' % (celltype, ii)
        model = moose.Neutral('/model/%s' % (test_id))
        data = moose.Neutral('/data/%s' % (test_id))    
        cell = cellclass(cellclass.prototype, '%s/%s' % (model.path, celltype))
        stimulus = moose.PulseGen('%s/stimulus' % (model.path))
        stimulus.firstDelay = DELAY
        stimulus.firstWidth = WIDTH
        stimulus.firstLevel = low + ii * step
        stimulus.connect('outputSrc', cell.soma, 'injectMsg')
        vm_table = moose.Table('%s/Vm' % (data.path))
        vm_table.stepMode = moose.TAB_BUF
        vm_table.connect('inputRequest', cell.soma, 'Vm')
        spike_table = moose.Table('%s/spiketimes' % (data.path))
        spike_table.connect('inputRequest', cell.comp[cell.presyn], 'Vm')
        stim_table = moose.Table('%s/stimulus' % (data.path))
        stim_table.stepMode = moose.TAB_BUF
        stim_table.connect('inputRequest', stimulus, 'output')
        rec_info[test_id] = {'Vm': vm_table,    
                             'spike': spike_table,
                             'stimulus': stim_table,
                             'current': stimulus.firstLevel}

    moose.context.setClock(0, simdt)
    moose.context.setClock(1, simdt)
    moose.context.setClock(2, simdt)
    moose.context.setClock(3, simdt)
    moose.context.setClock(4, simdt)
    moose.context.setClock(5, plotdt)
    moose.context.useClock(5, '/data/##[ISA=Table]')
    if config.clockjob.autoschedule != 1:
        moose.context.useClock(0, '/model/##[TYPE=Compartment]', 'init')
        moose.context.useClock(1, '/model/##[TYPE=Compartment]', 'process')
        moose.context.useClock(2, '/model/##[TYPE!=Compartment]')
    moose.context.reset()
    moose.context.step(SIMTIME)
    for test_id, info in rec_info.items():
        grp = h5file.require_group(test_id)
        grp.attrs['current'] = info['current']
        vm = grp.create_dataset('Vm', shape=(len(info['Vm']),), 
                                data=np.asarray(info['Vm']),
                                dtype=np.float64,
                                compression='gzip')
        vm.attrs['dt'] = plotdt
        stim = grp.create_dataset('stimulus', shape=(len(info['stimulus']),), 
                                data=np.asarray(info['stimulus']),
                                dtype=np.float64,
                                compression='gzip')
        stim.attrs['dt'] = plotdt 
        spike = grp.create_dataset('spike', shape=(len(info['spike']),),
                                   data=np.asarray(info['spike']),
                                   dtype=np.float64,
                                   compression='gzip')
                                  
I_START = -2000e-12
I_END = 2000e-12
I_STEP = 100e-12

if __name__ == '__main__':
    fname = 'f_i_tests_{}.h5'.format(sys.argv[1])
    with h5.File(fname) as h5file:
        print 'Running current step test for', sys.argv[1]
        h5file.attrs['TIMESTAMP'] = config.timestamp.isoformat()
        current_step_test(sys.argv[1], I_START, I_END, I_STEP, h5file)

    print 'Finished. Data in', fname


