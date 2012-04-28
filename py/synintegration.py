# synintegration.py --- 
# 
# Filename: synintegration.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Apr 24 15:30:05 2012 (+0530)
# Version: 
# Last-Updated: Sat Apr 28 18:41:44 2012 (+0530)
#           By: subha
#     Update #: 433
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This script is to set up a single spiny stellate cell and deliver
# artificial synaptic inputs to study the synaptic integration time. I
# noticed in network model that the spinystellate cells were
# responding about 40 ms after a stimulus. The purpose of this
# simulation is to find out the reason for that.
# 
# 

# Code:
import sys

sys.path.append('/data/subha/chamcham_moose/python')
import numpy as np
import h5py as h5
import moose
from trbsim import Simulation
from spinystellate import SpinyStellate
from tcr import TCR
from deepbasket import DeepBasket
from trbnetdata import TraubFullNetData
from datetime import datetime

netdata = TraubFullNetData()

def get_syninfo(filename, cellname, srctype):
    """Load the synapses on cell from network file specified by
    filename"""
    ret = None
    with h5.File(filename, 'r') as netfile:
        syntab = netfile['/network/synapse'][:]
        indices = np.nonzero(np.char.startswith(syntab['dest'], cellname) & np.char.startswith(syntab['source'], srctype))[0]
        subarray = syntab[indices]
        print subarray
        comps = [int(entry[-1]) for entry in np.char.rpartition(subarray['dest'], '_')]
        print comps
        gbar = subarray['Gbar'].astype(float)
        tau1 = subarray['tau1'].astype(float)
        tau2 = subarray['tau2'].astype(float)
        Ek =  subarray['Ek'].astype(float)
        ret = (comps, subarray['type'], gbar, tau1, tau2, Ek)
    return ret

def synintegration2(filename, target_cell):
    raise NotImplementedError
    syntab = None
    netfile = h5.File(filename, 'r')
    syntab = netfile['network/synapse'][:]
    idx = np.char.startswith(syntab['dest'], target_cell)
    syntab = syntab[idx]
    sim = Simulation('synintegration')
    cell = SpinyStellate(SpinyStellate.prototype, 
                         sim.model.path + '/SpinyStellate')
    vm_table = cell.comp[cell.presyn].insertRecorder('Vm', 'Vm', sim.data)
    # Create a common spike gen object
    sp = moose.SpikeGen('spike', sim.model)
    sp.threshold = 0.0
    sp.edgeTriggered = 1
    sptab = moose.Table('spike', sim.data)
    sptab.stepMode = 3
    sptab.connect('inputRequest', sp, 'state')
    for ii in range(len(syntab)):
        source = syntab['source'][ii].partition('/')[0]
        dest, compname = syntab['dest'][ii].split('/')
        if source == dest:
            print 'source = dest', source, filename
            continue
        comp_no = compname.split('_')[-1]
        comp = cell[int(comp_no)]
        weight = 1.0
        if syntab['type'][ii] == 'nmda':
            chan = moose.NMDAChan('nmda_from_%s' % 
                                  (source.split('_')[0]), comp)
            chan.MgConc = 1.5            
            weight = syntab['Gbar'][ii]
        else:
            chan = moose.SynChan('%s_from_%s' % 
                                 (syntab['type'][ii], source.split('_')[0]), 
                                 comp)
        chan.Gbar = syntab['Gbar'][ii]
        chan.tau1 = syntab['tau1'][ii]
        chan.tau2 = syntab['tau2'][ii]
        chan.Ek = syntab['Ek'][ii]
        comp.connect('channel', chan, 'channel')
        sp.connect('event', chan, 'synapse')
        count = chan.numSynapses
        if source_type == 'TCR':
            chan.delay[count-1] = 1e-3 # thalamocortical delay
        else:
            chan.delay[count-1] = 0.05e-3 # default delay
        chan.weight[count-1] = weight
        gktable = moose.Table('%s_%s_%s' % (cell.name, comp.name, chan.name), sim.data)
        gktable.stepMode = 3
        gktable.connect('inputRequest', chan, 'Gk')
            

thalamic = ['nRT', 'TCR']

def delay(pretype, posttype):
    if pretype in thalamic and posttype not in thalamic:
        return 1e-3
    elif pretype not in thalamic and posttype in thalamic:
        return 5e-3
    else:
        return 0.05e-3
    
def ampa(pre_ix, post_ix):
    """Return parameters for AMPA synapse."""
    tau = netdata.tau_ampa[pre_ix][post_ix]
    return {'tau1': tau,
            'tau2': tau,
            'Gbar': netdata.g_ampa_baseline[pre_ix][post_ix] * \
                tau * 1e3 / np.e,
            'Ek': 0.0,
            'weight': 1.0,
            'classname': 'SynChan'}
    
def gaba(pre_ix, post_ix):
    """Return parameters for GABA synapse with a single component."""
    return {'tau1': netdata.tau_gaba[pre_ix][post_ix],
            'tau2': 0.0,
            'Gbar': netdata.g_gaba_baseline[pre_ix][post_ix],
            'Ek': -81e-3,
            'weight': 1.0,
            'classname': 'SynChan'}

def nmda(pre_ix, post_ix):
    """Return parameters for NMDA synapse."""
    return {'tau1': netdata.tau_nmda[pre_ix][post_ix],
            'tau2': 0.0,
            'Gbar': netdata.g_nmda_baseline[pre_ix][post_ix],
            'weight': netdata.g_nmda_baseline[pre_ix][post_ix],
            'Ek': 0.0,
            'classname': 'NMDAChan'}


def make_synapse(precell, postcell, syntype):
    pretype = precell.__class__.__name__
    posttype = postcell.__class__.__name__
    pre_ix = netdata.celltype.index(pretype)
    post_ix = netdata.celltype.index(posttype)
    precomp = precell.comp[precell.presyn]
    try:
        postcomp_ix = netdata.allowed_comps[pre_ix][post_ix][0]
        postcomp = postcell.comp[postcomp_ix]
    except IndexError:
        return None
    kwargs = {}
    classname = 'SynChan'
    if syntype == 'ampa':
        kwargs = ampa(pre_ix, post_ix)
    elif syntype == 'nmda':
        kwargs = nmda(pre_ix, post_ix)
    else:
        kwargs = gaba(pre_ix, post_ix)
    kwargs['delay'] = delay(netdata.celltype[pre_ix], 
                            netdata.celltype[post_ix])
    kwargs['name'] = '%s_from_%s' % (syntype, pretype)    
    kwargs['classname'] = classname
    if kwargs['Gbar'] == 0.0:
        return None
    print '== Synapse:', precell.path, '<=', postcell.path, '=='
    print kwargs
    print '-- !Synapse Args --'
    return precomp.makeSynapse(postcomp, **kwargs)

def record_data(datacontainer, tabname, targetobj, targetfield):
    tab = moose.Table(tabname, datacontainer)
    tab.stepMode = 3
    tab.connect('inputRequest', targetobj, targetfield)
    return tab
    
def threecell_test():
    sim = Simulation('threecell')
    tcr_ix = netdata.celltype.index('TCR')
    ss_ix = netdata.celltype.index('SpinyStellate')
    basket_ix = netdata.celltype.index('DeepBasket')
    ss = SpinyStellate(SpinyStellate.prototype,
                       sim.model.path + '/SpinyStellate')
    tcr = TCR(TCR.prototype, sim.model.path + '/TCR')
    basket = DeepBasket(DeepBasket.prototype, sim.model.path + '/DeepBasket')
    cells = [ss, tcr, basket]
    for postcell in cells:
        for precell in cells:
            if precell == postcell:
                continue
            for syntype in ['ampa', 'nmda', 'gaba']:
                syn = make_synapse(precell, postcell, syntype)
                if syn is not None:
                    tab = record_data(sim.data, 'gk_%s_%s' % (postcell.name, syn.name), syn, 'Gk')
        record_data(sim.data, 'Vm_%s' % (postcell.name), postcell.comp[postcell.presyn], 'Vm')
    pulsegen = moose.PulseGen('pulse', sim.model)
    pulsegen.firstDelay = 3.0
    pulsegen.firstLevel = 1e-9
    pulsegen.firstWidth = 2e-3
    pulsegen.trigMode = moose.FREE_RUN
    pulsegen.connect('outputSrc', tcr.soma, 'injectMsg')
    ptable = moose.Table('pulse', sim.data)
    ptable.stepMode = 3
    ptable.connect('inputRequest', pulsegen, 'output')
    sim.schedule(simdt=25e-6, plotdt=25e-5)
    for t in range(10):
        start = datetime.now()
        sim.run(1.0)
        end = datetime.now()
        delta = end - start
        print 'Ran', t+1, '-th second in', delta.seconds + 1e-6 * delta.microseconds
    sim.save_data_h5('threecells.h5')
    
    

def test_synintegration(filename, target_cell, source_type):
    sim = Simulation('synintegration')
    cell = SpinyStellate(SpinyStellate.prototype, 
                         sim.model.path + '/SpinyStellate')
    vm_table = cell.comp[cell.presyn].insertRecorder('Vm', 'Vm', sim.data)
    # Create a common spike gen object
    sp = moose.SpikeGen('spike', sim.model)
    sp.threshold = 0.0
    sp.edgeTriggered = 1
    sptab = moose.Table('spike', sim.data)
    sptab.stepMode = 3
    sptab.connect('inputRequest', sp, 'state')
    (comp_indices, syntype, gbar, tau1, tau2, Ek, ) = \
        get_syninfo(filename, target_cell, source_type)    
    # Set up the synapses
    for ii in range(len(comp_indices)):
        print '%d\t%s\t%g\t%g\t%g\t%g' % (comp_indices[ii], 
                                          syntype[ii], 
                                          gbar[ii], 
                                          tau1[ii], 
                                          tau2[ii], 
                                          Ek[ii])
        comp = cell.comp[comp_indices[ii]]
        weight = 1.0
        if syntype[ii] == 'nmda':
            chan = moose.NMDAChan('nmda_from_%s' % (source_type), comp)
            chan.MgConc = 1.5            
            weight = gbar[ii]
        else:
            chan = moose.SynChan('%s_from_%s' % (syntype[ii], source_type), 
                                 comp)
        chan.Gbar = gbar[ii]
        chan.tau1 = tau1[ii]
        chan.tau2 = tau2[ii]
        chan.Ek = Ek[ii]
        comp.connect('channel', chan, 'channel')
        sp.connect('event', chan, 'synapse')
        count = chan.numSynapses
        if source_type == 'TCR':
            chan.delay[count-1] = 1e-3 # thalamocortical delay
        else:
            chan.delay[count-1] = 0.05e-3 # default delay
        chan.weight[count-1] = weight
        gktable = moose.Table('%s_%s_%s' % (cell.name, comp.name, chan.name), sim.data)
        gktable.stepMode = 3
        gktable.connect('inputRequest', chan, 'Gk')
    pulsegen = moose.PulseGen('pulse', sim.model)
    pulsegen.firstDelay = 3.0
    pulsegen.firstLevel = 1.0
    pulsegen.firstWidth = 1e-3
    pulsegen.trigMode = moose.FREE_RUN
    pulsegen.connect('outputSrc', sp, 'Vm')
    ptable = moose.Table('pulse', sim.data)
    ptable.stepMode = 3
    ptable.connect('inputRequest', pulsegen, 'output')
    sim.schedule(simdt=1e-6, plotdt=1e-6)
    sim.run(10.0)
    sim.save_data_h5(filename.replace('network_', 'synintegration_'))

if __name__ == '__main__':
    # test_synintegration('/data/subha/rsync_ghevar_cortical_data_clone/2012_01_18/network_20120118_142820_7865.h5.new',
    #                     'SpinyStellate_0',
    #                     'TCR')
    threecell_test()
            
# 
# synintegration.py ends here
