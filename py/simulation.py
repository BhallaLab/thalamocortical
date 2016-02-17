# simulation.py --- 
# 
# Filename: simulation.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Thu Apr 30 02:25:06 2009 (+0530)
# Version: 
# Last-Updated: Tue Jan  1 09:33:57 2013 (+0530)
#           By: subha
#     Update #: 482
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# A manager class to handle a whole simulation.
# 
# 

# Change log:
#
# 2009-04-30 02:25:59 (+0530) - factored out the Simulation class from
# test.py to a separate simulation.py.
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
import os
from datetime import datetime
import tables
import numpy
import moose
import pymoose
import config
import smtplib

channel_lib = {}


class Simulation:
    """This class is a wrapper to control a whole simulation."""
    def __init__(self, name):
        self.name = name
        self.model = moose.Neutral('model')
        self.data = moose.Neutral('data')
        self.start_t = None
        self.end_t = None
        self.simdt = config.simdt
        self.plotdt = config.plotdt
        self.gldt = config.gldt
        self.simtime = 50e-3

    def schedule(self, simdt=None, plotdt=None, gldt=None):
        """Schedule the objects and do a reset.

        Fine control over scheduling of instances differenct classes.
        
        """
        config.LOGGER.info('Starting scheduling and reset')
        if simdt:
            self.simdt = simdt
        if plotdt:
            self.plotdt = plotdt
        if gldt:
            self.gldt = gldt
        config.context.setClock(0, self.simdt)
        config.context.setClock(1, self.simdt)
        config.context.setClock(2, self.simdt)
        config.context.setClock(3, self.simdt)
        config.context.setClock(4, self.simdt) # This is special clock for recording stimulus, which has a pulsewidth less than plotdt
        config.context.setClock(5, self.plotdt)
        
        if config.clockjob.autoschedule is 1:
            config.BENCHMARK_LOGGER.info('Automatic scheduling...')
            start = datetime.now()
            config.context.useClock(4, self.data.path + '/##[TYPE=Table]')
            config.context.useClock(5, self.data.path + '/stimulus/#[TYPE=Table]')
            end = datetime.now()
            delta = end - start
            config.BENCHMARK_LOGGER.info('Time to schedule Table objects: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
            start = datetime.now()
            config.context.reset()
            end = datetime.now()
            delta = end - start
            print '%%%%%%%%%%%%%%% reset'
            config.BENCHMARK_LOGGER.info('reset time: %g s' % (delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))
            # config.LOGGER.debug('*** Heap usage  after reset ***')
            # config.LOGGER.debug(config.heapy.heap())        
            return
        
        # start = datetime.now()
        # config.context.useClock(0, self.model.path + '/##')
        # end = datetime.now()
        # delta = end - start
        # config.BENCHMARK_LOGGER.info('Time to schedule all objects on default clock 0: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        # config.LOGGER.debug('*** Heap usage  beafore scheduling clock#0 ***')
        # config.LOGGER.debug(config.heapy.heap()) 
        start = datetime.now()
        # clock0_paths = ['%s/##[TYPE=%s' % (self.model.path, typename) for typename in ['Compartment', 'CaConc', 'HHChannel', 'HHGate', 'SynChan', 'NMDAChan', 'SpikeGen', 'RandomSpike', 'PulseGen']]
        # config.context.useClock(0, '%s/##[TYPE=Compartment],%s/##[TYPE=CaConc],%s/##[TYPE=HHCHannel],%s/##[TYPE=HHGate],%s/##[TYPE=SynChan]' % (self.model.path, self.model.path, self.model.path, self.model.path, self.model.path))
        config.context.useClock(0, '%s/##[TYPE=Compartment]' % (self.model.path), 'init')
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Time to schedule all model elements: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        start = datetime.now()
        config.context.useClock(1, self.model.path + '/##[TYPE=Compartment]', 'process')
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Time to schedule init method of Compartment objects: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))

        start = datetime.now()
        config.context.useClock(3, self.model.path + '/##[TYPE!=Compartment]')
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Time to schedule channel, caconc, objects: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        start = datetime.now()
        config.context.useClock(4, self.data.path + '/##[TYPE=Table]')
        end = datetime.now()
        # config.LOGGER.debug('*** Heap usage  after scheduling clock#3 ***')
#        config.LOGGER.debug(config.heapy.heap())        
        delta = end - start
        config.BENCHMARK_LOGGER.info('Time to schedule Table objects: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        config.context.useClock(5, self.data.path + '/stimulus/#[TYPE=Table]')
        start = datetime.now()
        config.context.reset()
        end = datetime.now()
        # config.LOGGER.debug('*** Heap usage  after reset ***')
#        config.LOGGER.debug(config.heapy.heap())        
        delta = end - start
        config.BENCHMARK_LOGGER.info('reset time: %g s' % (delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))
        if config.solver == 'hsolve':
            # HSolve gets created only after reset - so we have to schedule it afterwards.
            start = datetime.now()
            config.context.useClock(2, self.model.path +'/##[TYPE=HSolve]')
            end = datetime.now()
            delta = end - start
            config.BENCHMARK_LOGGER.info('Time to schedule HSolve objects: %g s' % (delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))
        # config.LOGGER.debug('*** Heap usage  after scheduling hsolve objects on Clock#2 ***')
#        config.LOGGER.debug(config.heapy.heap())        

    def run(self, time=None):
        """Run simulation for given time."""
        if time is None or time <= 0.0:
            time = self.simtime
        else:
            self.simtime = time
        config.LOGGER.info('Starting simulation for %g s' % (self.simtime))
        self.start_t = datetime.now()
        config.context.step(float(time))
        self.end_t = datetime.now()
        delta = self.end_t - self.start_t
        config.BENCHMARK_LOGGER.info('%s : RUNTIME for %g s of simulation: %g s' %  (self.name, 
                                                                                     self.simtime, 
                                                                                     delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds))

    def reset_and_run(self, simtime=None, simdt=None, plotdt=None, gldt=None, updateinterval=0.5):
        """Reset and run together. This function uses simtime and
        plotdt together ensure prior allocation of table vector."""
        if (simtime is None) or (simtime <= 0.0):
            simtime = self.simtime
        else:
            self.simtime = simtime
        # The following explicit tablesize setting is to avoid excess
        # memory allocation for table vectors.
        xdivs = int(self.simtime / self.plotdt + 0.5) + 1
        for table_id in self.data.children('/##[TYPE=Table]'):
            table = moose.Table(table_id)
            if table.stepMode == moose.TAB_SPIKE:
                continue
            table.xmin = 0.0
            table.xmax = self.simtime
            table.xdivs = xdivs
        self.schedule(simdt=simdt, plotdt=plotdt, gldt=gldt)
        time = updateinterval
        while (time <= simtime):
            self.run(updateinterval)
            config.LOGGER.info('Simulation done till: %f s of %f' % (time, simtime))
            time += updateinterval
        rest = simtime - (time-updateinterval)
        if rest > 0:
            self.run(rest)
        config.LOGGER.info('Finished simulation for: %f s' % (simtime))

    def dump_data(self, directory, time_stamp=False):
        """Save the data in directory. The files are prefixed with
        config.timestamp in HH.MM.SS_ format if time_stamp is True."""
        path = directory
        tables = []
        if not os.access(path, os.W_OK):
            config.LOGGER.warning('data directory: ' + path + ' is not writable')
            return
        else:
            if not os.access(path, os.F_OK):
                os.mkdir(path)
        if time_stamp:
            prefix = config.timestamp.strftime('%H.%M.%S') + '_'
        else:
            prefix = ''
        for table_id in self.data.children():
            table = moose.Table(table_id)
                

            tables.append(table)
            file_name = prefix + table.name + '.plot'
            file_path = os.path.join(path, file_name)
            if table.stepMode == 3:
                ts = numpy.linspace(0, self.simtime, len(table))
                numpy.savetxt(file_path, zip(ts, numpy.asarray(table)))
            else:
                table.dumpFile(file_path)
            config.LOGGER.info('Dumped data in %s' % (file_path))
        return tables

    def save_data_h5(self, filename, notes=''):
        """Save all the data tables in an hdf5 file.

        The structure in the hdf5 file mirrors the structure in the
        model under /data element - but only two levels deep."""
        config.LOGGER.debug('Start saving the data')
        starttime =  datetime.now()
        compression_filter =  tables.Filters(complevel=9, complib='zlib', fletcher32=True)
        h5file =  tables.openFile(filename,  mode = 'w',  title = 'Traub Network: timestamp: %s' % (config.timestamp.strftime('%Y-%m-%d %H:%M:%S')),  filters = compression_filter)
        h5file.root._v_attrs.simtime = self.simtime
        h5file.root._v_attrs.simdt = self.simdt
        h5file.root._v_attrs.plotdt = self.plotdt
        h5file.root._v_attrs['notes'] = notes
        h5file.root._v_attrs.timestamp = config.timestamp.strftime('%Y-%m-%d %H:%M:%S')
        # Save simulation configuration data. I am saving it both in
        # data file as well as network file as often the data file is
        # too large and may not be available if the simulation is
        # cancelled midway.
        runconfig = h5file.createGroup(h5file.root, 'runconfig', 'Simulation settings')
        for section in config.runconfig.sections():
            table_data = config.runconfig.items(section)
            if table_data:
                sectiontab = h5file.createTable(runconfig, section, table_data)
        for child_id in self.data.children():
            child = moose.Neutral(child_id)
            if child.className == 'Neutral':
                group = h5file.createGroup(h5file.root, child.name)
                for tab_id in child.children():
                    tab = moose.Table(tab_id)
                    print 'Saving', tab.path, 'of length:', len(tab)
                    if tab.className == 'Table':
                        carray = h5file.createCArray(group, tab.name, tables.FloatAtom(), shape=(len(tab),))
                        carray[:] = numpy.array(tab)
            elif child.className == 'Table':
                tab = moose.Table(child_id)
                carray = h5file.createCArray(h5file.root, tab.name, tables.FloatAtom(), shape=(len(tab),))
                carray[:] = numpy.array(tab)
                
                # NOTE: I found that saving the first value with first
                # forwards diffs saves a lot of space on
                # compression. But will go with direct method now to
                # see if it is worth it.
            else:
                raise Warning('Element %s of unhandled type %s under data element:  -- not known for saving.' % (tab.name, tab.className))
        h5file.close()
        endtime = datetime.now()
        delta = endtime - starttime
        config.BENCHMARK_LOGGER.info('Saved data to %s: %g s' % (filename, delta.days * 86400 + delta.seconds + delta.microseconds * 1e-6))


def test_simulation():
    sim = Simulation('test_simulation')

if __name__ == '__main__':
    sim = Simulation('test')
    send_email('target@gmail.com',
                   'sender@gmail.com',
                   'password')

# 
# simulation.py ends here
