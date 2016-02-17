# trbsim.py --- 
# 
# Filename: trbsim.py
# Description: 
# Author: 
# Maintainer: 
# Created: Tue Nov 23 18:35:45 2010 (+0530)
# Version: 
# Last-Updated: Mon Jul  1 12:05:24 2013 (+0530)
#           By: subha
#     Update #: 751
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This is the main simulation setup and execution code. This is the
# driver of the simulation.
# 
# 

# Change log:
# 
# 2011-09-07 16:12:23 (+0530) -- Switching to ConfigParser for finer
# control over the simulation. The configuration information will be
# copied over to the network and data files.
# 
# 

# Code:


import sys
from getopt import getopt
from subprocess import call
from collections import defaultdict
import numpy
from operator import itemgetter
from itertools import chain

import moose
import config

from simulation import Simulation
from trbnet import *
from trbutil import send_email
from kchans import KChannel

def get_targets(cell_paths):
    """Get the cells post synaptic to cells listed in cell_paths."""
    # print '###', cell_paths
    target_sources_map = defaultdict(list)
    for cp in cell_paths:
        # print '###', cp
        sg_list = moose.context.getWildcardList(cp + '/##[TYPE=SpikeGen]', True)
        for sg in sg_list:
            sg = moose.SpikeGen(sg)
            for synId in sg.neighbours('event', 2):
                print synId.path()
                syn = moose.Neutral(synId)
                if not(syn.className == 'SynChan' or syn.className == 'NMDAChan'):
                    continue
                for comp in moose.Neutral(synId).neighbours('channel', 2):
                    # print '-->', comp.path()
                    target_sources_map[comp.path().rpartition('/')[0]].append(cp)
    return target_sources_map

def get_max_targets(target_sources_map, frac=0.1):
    celltype_dict = defaultdict(dict)
    for target, sources in target_sources_map.items():
        celltype_dict[target.rpartition('/')[-1].rpartition('_')[0]][target] = len(sources)
    ret = defaultdict(list)
    for celltype, tsdict in celltype_dict.items():
        sortedcells = sorted(tsdict.items(), key=itemgetter(1))
        ll = int(frac * len(sortedcells))
        if ll == 0:
            ll = 1
        ret[celltype] = sortedcells[-ll:]
    return ret
        
def main(*argv):
    """-c cell-count-file

    -g celltype-graph-file

    -s synaptic-scaling-file

    -v synaptic-conductance-values-file

    -t simulationtime

    -d simdt

    -p plotdt

    -r randomseed

    -s -- do explicit scheduling

    --reseed -- seed the rng with PID of current process

    -x prelease_unknown -- default release probability in stochastic simulations.

    --stochastic -- synapses are stochastic

    -i -- interstimulus interval
    """
    print '******** Starting run _{timestamp}_{PID} :', config.filename_suffix
    config.LOGGER.info('Running moose version: %s revision: %s' % (moose.version, moose.revision))
    celltypegraph_file = None
    cellcount_file = None
    synscale_file = None
    synval_file = None
    notes = None
    netfile = None
    optlist, args = getopt(argv[0], 'i:m:n:ehc:s:v:t:d:p:r:g:lx:', ['reseed', 'stochastic', 'netfile='])
    print 'Program arguments: '
    print 'Options:', optlist
    print 'Args:', args
    for arg in optlist:        
        if arg[0] == '-c':
            cellcount_file = arg[1]
        elif arg[0] == '-n':
            notes = arg[1]
        elif arg[0] == '-m':
            with open(arg[1], 'r') as notes_file:
                notes = '\n'.join(notes_file.readlines())
        elif arg[0] == '-s':
            synscale_file = arg[1]
        elif arg[0] == '-v':
            synval_file = arg[1]
        elif arg[0] == '-t':
            config.simtime = float(arg[1])
        elif arg[0] == '-d':
            config.simdt = float(arg[1])
        elif arg[0] == '-p':
            config.plotdt = float(arg[1])
        elif arg[0] == '-r':
            config.rngseed = int(arg[1])
        elif arg[0] == '-g':
            celltypegraph_file = arg[1]
        elif '--reseed' == arg[0]:
            config.to_reseed = True
        elif '--stochastic' == arg[0]:
            config.stochastic = True
        elif '--netfile' == arg[0]:
            netfile = arg[1]
        elif '-e' == arg[0]:
            config.solver = 'ee'
        elif arg[0] == '-l':
            config.clockjob.autoschedule = 0
        elif arg[0] == '-x':
            config.default_releasep = float(arg[1])
            assert((config.default_releasep < 1.0) and (config.default_releasep >= 0.0))
        elif arg[0] == '-i':
            config.runconfig.set('stimulus', 'isi', arg[1])
            print 'Set isi to', arg[1]
        elif '-h' == arg[0]:
            print main.__doc__
            sys.exit(0)

        
    sim = Simulation('sim%s' % (config.filename_suffix))
    net = TraubNet(container=sim.model, netfile=netfile)
    if notes is not None:
        net.tweaks_doc.append(notes)
    if celltypegraph_file is None:
        net._generate_celltype_graph()
    else:
        net.setup_from_celltype_file(celltype_file=celltypegraph_file)
    net.set_populations()
    net.set_conductances(synval_file)
    net.tune_conductances(synscale_file)
    scale_nmda = config.runconfig.get('NMDA', 'conductance_scale')
    if scale_nmda:
        net.scale_conductance('gnmda', {('*', '*'):float(scale_nmda)})
    scale_ampa = config.runconfig.get('AMPA', 'conductance_scale')
    if scale_ampa:
        net.scale_conductance('gampa', {('*', '*'):float(scale_ampa)})
    scale_gaba = config.runconfig.get('GABA', 'conductance_scale')
    if scale_gaba:
        net.scale_conductance('ggaba', {('*', '*'):float(scale_gaba)})    
    net.set_unknown_prelease(config.default_releasep)
    net._generate_cell_graph()
    net.create_network()
    if config.runconfig.get('stimulus', 'ectopic_spike') in ['1', 'yes', 'Yes', 'YES', 'true', 'True', 'TRUE']:
        net.setup_ectopic_input()
    # Comment out the following line to use original K+ reversal potential.
    # net.tweak_Ek(KChannel, 10e-3)
    net.randomize_passive_properties()
    net.randomize_active_conductances()
    background_celltype = config.runconfig.get('stimulus', 'background')
    probe_celltype = config.runconfig.get('stimulus', 'probe')
    bg_count = int(config.runconfig.get('stimulus', 'bg_count'))
    probe_count = int(config.runconfig.get('stimulus', 'probe_count'))
    onset = float(config.runconfig.get('stimulus', 'onset'))
    bg_interval = float(config.runconfig.get('stimulus', 'bg_interval'))
    bg_interval_spread = float(config.runconfig.get('stimulus', 'bg_interval_spread'))
    pulse_width = float(config.runconfig.get('stimulus', 'pulse_width'))
    isi = float(config.runconfig.get('stimulus', 'isi'))
    print 'trbsim - isi', isi
    amp = float(config.runconfig.get('stimulus', 'amplitude'))
    print bg_interval, pulse_width, isi
    if bg_interval_spread > 0.0:
        num_bg_pulses = (config.simtime - onset) / (bg_interval + bg_interval_spread/2)
    else:
        num_bg_pulses = 0
    stim_dict = net.setup_stimulus(bg_cells=background_celltype,
                       probe_cells=probe_celltype,
                       bg_count=bg_count,
                       probe_count=probe_count,
                       stim_onset=onset,
                       bg_interval=bg_interval,
                       bg_interval_spread=bg_interval_spread,
                       num_bg_pulses=num_bg_pulses,
                       pulse_width=pulse_width,
                       isi=isi,
                       level=amp)
    # Record post synaptic Vm for cells with highest convergent connections
    pspfrac = float(config.runconfig.get('record', 'pspfrac'))
    pspcontainer = moose.Neutral('%s/psp' % (sim.data.path))
    if pspfrac > 0:
        bg_target_cells = [comppath.rpartition('/')[0] for comppath in stim_dict['bg_targets']]
        target_source_map = get_targets(bg_target_cells)
        max_targets = get_max_targets(target_source_map, pspfrac)
        for celltype, cells in max_targets.items():
            config.LOGGER.debug('Recording Vm for 1 hop from stimulus for %s: %s' % (celltype, cells))
            for cellpath, count in cells:
                soma = MyCompartment(cellpath + '/comp_1')
                if soma.className == 'Compartment':
                    soma.insertRecorder(cellpath.rpartition('/')[-1], 'Vm', pspcontainer)
    # for cpath in stim_dict['bg_targets']:
    #     cell = cpath.rpartition('/')[0]
    #     comp = MyCompartment(cpath)
    #     comp.insertRecorder(
    #     , stim_dict['probe_targtes']):
        
    if net.from_netfile is not None and synscale_file is not None:
        net.scale_synapses(synscale_file)
    net.save_network_model(config.MODEL_FILENAME)
    net.verify_saved_model(config.MODEL_FILENAME)
    net.save_cell_network(config.MODEL_FILENAME + '.new')
    net.setup_spike_recording(sim.data)
    # net.setup_current_injection_test([-0.9e-9, -0.3e-9, 0.0, 0.1e-9, 0.0, 0.3e-9, 0.0, 0.9e-9], 1.0, 1.0, sim.data)
    # Setup the bias currents
    bias_current_delays = defaultdict(dict)
    bias_current_levels = defaultdict(dict)
    bias_current_widths = defaultdict(dict)
    for key, value in config.runconfig.items('bias_current'):
        if not value:
            continue
        tokens = key.split('_')
        celltype = tokens[0]
        param = tokens[1]
        index = int(tokens[2])
        try:
            if param == 'delay':
                bias_current_delays[celltype][index] = float(value)
            elif param == 'width':
                bias_current_widths[celltype][index] = float(value)
            elif param == 'level':
                bias_current_levels[celltype][index] = float(value)
        except ValueError, e:
            print 'Value:"', value, '"'
            raise(e)
    for key, level_dict in bias_current_levels.items():
        delay_dict = bias_current_delays[key]
        width_dict = bias_current_widths[key]
        level_list = []
        delay_list = []
        width_list = []
        for ii in range(len(level_dict)):
            level_list.append(level_dict[ii])
            delay_list.append(delay_dict[ii])
            width_list.append(width_dict[ii])
        net.setup_bias_current(key, level_list, delay_list, width_list, sim.data)
        
    for celltype, cellcount in config.runconfig.items('vm_record'):
        if cellcount is None:
            cellcount=0
        elif '.' in cellcount or 'e' in cellcount:
            cellcount = int(config.runconfig.get('cellcount', celltype)) * float(cellcount)
        net.setup_Vm_recording(sim.data, celltype, numcells=int(cellcount))
    if config.runconfig.get('record', 'lfp') in ['Yes', 'YES', 'yes', 'true', 'True', 'TRUE', '1']:
        electrode_depths = numpy.arange(200e-6, 2200e-6, 200e-6)
        for depth in electrode_depths:
            net.setup_lfp_recording('electrode_%dum' % (int(depth*1e6)), depth, sim.data) # lfp electrode at 200 um depth interval
        # We also want to separate cell class specific LFP.
        lfp_cellclasses = ['SupPyrRS', 'SupPyrFRB', 'TuftedIB', 'TuftedRS', 'NontuftedRS']
        for cellclass in lfp_cellclasses:
            net.setup_lfp_recording('electrode_%dum_%s' % (1000, cellclass), 1000e-6, sim.data, [cellclass]) 
            net.setup_lfp_recording('electrode_%dum_%s' % (2000, cellclass), 2000e-6, sim.data, [cellclass]) 
    ###############################
    # for debugging only
    #
    print 'PRINTING Synapse info'
    config.LOGGER.debug('PRINTING Synapse info')
    for celltype, count in config.runconfig.items('cellcount'):
        if int(count) == 0:
            continue
        synchans = moose.context.getWildcardList('%s/%s_0/##[ISA=SynChan]' % (net.network_container.path, celltype), True)
        nmdachans = moose.context.getWildcardList('%s/%s_0/##[TYPE=NMDAChan]' % (net.network_container.path, celltype), True)
        config.LOGGER.debug('No. of synchans on SpinyStellate_0 : %d' % (len(synchans)))
        config.LOGGER.debug('No. of nmdachans on SpinyStellate_0 : %d' % (len(nmdachans)))
        syn_tabs = []
        for synid in chain(synchans, nmdachans):
            synapse = moose.SynChan(synid)
            config.LOGGER.debug('Recording from %s' % (synapse.path))
            comp = moose.Neutral(synapse.parent)
            cell = moose.Neutral(comp.parent)
            synapse_data = moose.Neutral('synapse', sim.data)
            tab = moose.Table('gk_%s_%s_%s' % (cell.name, comp.name, synapse.name), synapse_data)
            tab.stepMode = 3
            config.LOGGER.info('Connected recording table %s: %s' % (tab.name, tab.connect('inputRequest', synapse, 'Gk')))
            syn_tabs.append(tab)
            tab = moose.Table('Ik_%s_%s_%s' % (cell.name, comp.name, synapse.name), synapse_data)
            tab.stepMode = 3
            config.LOGGER.info('Connected recording table %s: %s' % (tab.name, tab.connect('inputRequest', synapse, 'Ik')))
            syn_tabs.append(tab)
            config.LOGGER.debug('== start synapse info')
            config.LOGGER.debug(synapse.path)
            msg_list = pymoose.listmsg(synapse)
            for msg in msg_list:
                config.LOGGER.debug(msg)
            fields = pymoose.getfields(synapse)
            for field in fields.items():
                config.LOGGER.debug(field)
            for ii in range(synapse.numSynapses):
                config.LOGGER.debug('delay[%d]: %g' % (ii, synapse.delay[ii]))
                config.LOGGER.debug('weight[%d]: %g' % (ii, synapse.weight[ii]))
            config.LOGGER.debug('-- end synapse info')
        
    # debugging till here
    ###############################
    sim.reset_and_run(simtime=config.simtime, simdt=config.simdt, plotdt=config.plotdt)
    sim.save_data_h5(config.DATA_FILENAME, '\n'.join(net.tweaks_doc))
    # print '******** Finished run _{timestamp}_{PID} :', config.filename_suffix
     # Now we try to send an email alert using a gmail account
#    subject = 'Simulation started at: %s is over' % (config.timestamp.strftime('%Y-%m-%d %H:%M:%S'))
    # subject = 'Simulation is over. Data in' + config.DATA_FILENAME
    # try:
    #     print 'Sending email'
    #     send_email('ray.subhasis@gmail.com', 
    #                'subha.gspace@gmail.com',
    #                'subha2google',
    #                subject=subject,
    #                body=notes)
    # except Exception, e:
    #     print e
    try:
        call(['gzip',config.LOG_FILENAME]) # To save space compress the log file.
    except Exception, e:
        print "Error: gzip", config.LOG_FILENAME, 'failed. The logfile was not compressed.'
        print e
       
    
import pymoose
def op(node):
    if not isinstance(node, moose.Neutral):
        node = moose.Neutral(node)
    if node.className == 'PulseGen':
        print 'Found PulseGen object %s' % (node.path)
    if not node.className == 'Table':
        pymoose.showmsg(node)

def check_isi(start, stop):
    raise NotImplementedError('This will not work well')
    last_isi = stop
    current_isi = start
    while current_isi < last_isi:
        config.runconfig.set('stimulus', 'isi', str(current_isi))
        sim = TraubSim()
        sim.runsim()
        
import traceback
if __name__ == '__main__':
    subject = 'Simulation started at %s is over.' % (config.timestamp.strftime('%Y-%m-%d %H:%M:%S'))
    body = 'Finished successfully. Data in %s' % (config.DATA_FILENAME)
    try:
        main(sys.argv[1:])
    except Exception, e:
        body = '<br/>'.join(traceback.format_exc().split('\n'))
        traceback.print_exc(file=sys.stdout)
    finally:
        send_email('target@gmail.com',
                   'sender@gmail.com',
                   'password',
                   subject=subject,
                   body=body)

    # pymoose.df_traverse(moose.Neutral('/'), op)


# 
# trbsim.py ends here
