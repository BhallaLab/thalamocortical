# tcr_spikegen.py --- 
# 
# Filename: tcr_spikegen.py
# Description: 
# Author: 
# Maintainer: 
# Created: Mon Jul  1 11:43:20 2013 (+0530)
# Version: 
# Last-Updated: Tue Jul  2 16:11:57 2013 (+0530)
#           By: subha
#     Update #: 125
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



import moose
import utils
import config

simtime = 100e-3
simdt = 1e-6

class TCRSpikeGen(moose.SpikeGen):
    """Subclass of SpikeGen with the makeSynapse function like the
    specialized compartment in traub net."""
    def __init__(self, *args):
        moose.SpikeGen.__init__(self, *args)
        self.amplitude = 1.0
        self.edgeTriggered = True
        
    def makeSynapse(self, target,
                    classname='SynChan',
                    name='synapse',
                    threshold=0.0, 
                    absRefract=0.0, 
                    Ek=0.0, 
                    Gbar=None, 
                    tau1=None, 
                    tau2=None,
                    weight=1.0,
                    delay=0.0,
                    Pr=1.0):
        """This replicates the makeSynapse in MyCompartment class in
        compartment.py"""
        classobj = eval('moose.' + classname)
        synapse = classobj(name, target)
        synapse.Ek = float(Ek) # TODO set value according to original model
        synapse.Gbar = float(Gbar) # TODO set value according to original model
        synapse.tau1 = float(tau1)
        synapse.tau2 = float(tau2)
        target.connect('channel', synapse, 'channel')
        if not self.connect('event', synapse, 'synapse'):
            raise Exception('Error creating connection: %s->%s' % (spikegen.path, synapse.path))
        num_synapses = synapse.numSynapses
        synapse.delay[num_synapses - 1] = float(delay)
        synapse.weight[num_synapses - 1] = float(weight)
        if config.stochastic:
            synapse.initPr[num_synapses - 1] = Pr
        config.LOGGER.debug('Created synapse: %s of type %s' % (synapse.path, synapse.className))
        return synapse

    def insertRecorder(self, objname, fieldname, data_container):
        """Replaces insertRecorder in MyCompartment. Always connects
        'event' message to 'input' of Table"""
        if config.context.exists('%s/%s' % (data_container.path, objname)):
            raise Exception('Trying to recreate existing data recorder: %s/%s' % (data_container.path, objname))
        table = moose.Table(objname, data_container)
        table.stepMode = moose.TAB_SPIKE
        table.stepSize = 0.1
        self.connect('state', table, 'inputRequest')
        config.LOGGER.debug('connected %s.event to %s.input. Originally requested %s' % (self.path, table.path, fieldname))
        return table

    def insertCaRecorder(self, *args, **kwargs):
        config.LOGGER('Ca recording requested for %s - do nothing' % (self.__class__.__name__))
        return None
        


class TCR(moose.Cell):
    """This class replaces the TCR cell with TimeTable object"""
    presyn = 1 # this is the only element to be passed as a presynaptic component
    prototype = None
    num_comp = 1
    def __init__(self, *args):
        moose.Cell.__init__(self, *args)
        self.comp = [None, TCRSpikeGen('comp_1', self)]
        self.soma = self.comp[1] # This is to fool some functions which expect a soma to be in the cell
        config.LOGGER.info('Created spikegen %s' % (self.soma.path))

TCR.prototype = TCR('/library/TCR')


if __name__ == '__main__':
    import pylab

    model = moose.Neutral('/model')
    data = moose.Neutral('/data')
    cell = TCR(TCR.prototype, 'TCR_0', model)
    vm = cell.soma.insertRecorder('TCR_0_Vm', 'Vm', data)
    comp = moose.Compartment('/model/comp')
    comp.Em = -70e-3
    comp.initVm = -70e-3
    comp.Rm = 1e9
    comp.Cm = 1e-12
    comp_vm = moose.Table('/data/comp_Vm')
    comp_vm.stepMode = 3
    comp_vm.connect('inputRequest', comp, 'Vm')
    cell.soma.makeSynapse(comp, Gbar=1e-9, tau1=5e-3, tau2=130e-3, delay=0.1e-3)
    gk_table = moose.Table('/data/Gk')
    gk_table.stepMode = 3
    gk_table.connect('inputRequest', moose.SynChan(comp.path + '/synapse'), 'Gk')
    pulsegen = moose.PulseGen('/model/pulse')
    pulsegen.firstDelay = 10e-3
    pulsegen.firstWidth = 20e-3
    pulsegen.firstLevel = 1e-9
    pulsegen.connect('outputSrc', cell.soma, 'Vm')
    pulse_table = moose.Table('/data/pulse')
    pulse_table.stepMode = 3
    pulsegen.connect('outputSrc', pulse_table, 'input')
    cell.soma.threshold = 0.5e-9
    config.context.setClock(0, simdt, 0)
    config.context.useClock(0, '%s/##' % (model.path))
    config.context.useClock(0, '%s/##' % (data.path))
    config.context.reset()
    config.context.step(simtime)    
    vm_arr = pylab.array(vm)
    # vm_arr = pylab.unique(vm_arr)
    pylab.plot(vm_arr, pylab.ones(len(vm_arr)), 'x', label='spike event')
    ts = pylab.linspace(0, simtime, len(pylab.array(comp_vm)))
    pylab.plot(ts, pylab.array(comp_vm)*1e3, label='post synaptic Vm (mV)')
    pylab.plot(ts, pylab.array(pulse_table)*1e9, label='stimulus (nA)')
    pylab.plot(ts, pylab.array(gk_table)*1e9, label='Gk (nS)')
    pylab.legend()
    pylab.show()
    

# 
# tcr_spikegen.py ends here
