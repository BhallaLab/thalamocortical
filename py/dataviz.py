#!/usr/bin/env python

# Filename: dataviz.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Wed Dec 15 10:16:41 2010 (+0530)
# Version: 
# Last-Updated: Mon Feb 14 16:52:38 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 636
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This is for visualizing neuronal activity in animation from a hdf5
# data file.
# 
# Decided to use matplotlib/mlab instead of mayavi for the sake of ease of coding.

# Change log:
# 
# 2010-12-15 10:17:49 (+0530) -- initial version
#
# 2010-12-17 11:30:12 (+0530) -- working matplotlib 2D animation with
# randomly generated numbers.
#
# 2010-12-21 11:53:32 (+0530) -- realized that a better way to
# organize data would be to create /data/spike /data/Vm and /data/Ca
# in the MOOSE model and the corresponding tables under those with
# same name as the cell it is recording from. Depending on table name
# suffix is as bad as filename extensions in Windows - one has to be
# consistent with the assumptions about table names between the
# simulation code and the data analysis code.  Hence forking this away
# into code for analyzing newer data.
#
# 2011-02-11 15:26:02 (+0530) -- scrapped matplotlib/mayavi approach
# and going for simple 2D rasters with option for selecting tables and
# scrolling (using Qt).
#

# Code:

import os
import sys
import time
import re

import numpy as np
from pysparse.sparse.spmatrix import ll_mat
import tables
from PyQt4 import QtCore, QtGui
from PyQt4 import Qwt5 as Qwt


ITERS = 1000

def train_to_cellname(trainname):
    """Extract the name of the generator cell from a spike/Vm/Ca train name"""
    return trainname[:trainname.rfind('_')]

class TraubData:
    """Data model for exploring simulation data.

    datafilename -- path of the data file

    h5f -- hdf5 file object
    """
    def __init__(self, filename):
        self.datafilename = filename
        self.h5f = None
        self.spiketrain_names = []
        self.Vm_names = []
        self.Ca_names = []
        self.simtime = None
        self.simdt = None
        self.plotdt = None
        self.timestamp = None
        
        try:
            self.h5f = tables.openFile(self.datafilename)            
        except Exception, e:
            raise e
        try:
            self.simtime = self.h5f.root._v_attrs.simtime
            self.simdt = self.h5f.root._v_attrs.simdt
            self.plotdt = self.h5f.root._v_attrs.plotdt
        except AttributeError, e:
            print e
        
        try:
            for train in self.h5f.root.spikes:
                self.spiketrain_names.append(train.name)
        except tables.NoSuchNodeError:
            print 'No node called /spikes'
        try:
            for vm_array in self.h5f.root.Vm:
                self.Vm_names.append(vm_array.name)
        except tables.NoSuchNodeError:
            print 'No node called /Vm'

        try:
            for ca_array in self.h5f.root.Ca:
                self.Ca_names.append(cell_name)
        except tables.NoSuchNodeError:
            print 'No node called /Ca'

            
    def get_data_by_name(self, names, datatype):
        """Return the data tables listed in names.

        If datatype is None then use the name to recognize type -
        ending with ca will be a [Ca2+] table, ending with Vm will be
        a Vm table, with spike will be a spike table.
        
        """
        ret = None
        target_node = None
        if datatype.lowercase() == 'ca':
            target_node = '/Ca'
        elif datatype.lowercase() == 'vm':
            target_node = '/Vm'
        elif datatype.lowercase() == 'spike':
            target_node = '/spikes'
        else:
            print 'Invalid data type: %s' % (datatype)
            return None
        ret = [self.h5f.getNode(target_node + '/' + name) for name in names]
        return ret

    def get_data_by_re(self, pattern, datatype):
        """Returns data for cells whose name match regular expression 'pattern'"""
        ret = []
        regex = re.compile(pattern)
        if datatype == 'spike':
            target_node = '/spikes'
        elif datatype == 'Vm' or datatype == 'vm':
            target_node = '/Vm'
        elif datatype == 'Ca' or datatype == 'ca':
            target_node = '/Ca'
        else:
            raise Exception('Invalid datatype specified: %s' %(datatype))
                            
        try:
            ret = [data for data in self.h5f.getNode(target_node) if regex.match(data.name)]
        except tables.NoSuchNodeError:
            print 'No such node: %s' % (target_node)
        return ret
                
    def __del__(self):
        # I think this is somewhat redundant - pytables takes care of
        # closing remaining open files.
        print 'del called'
        if self.h5f and self.h5f.isopen:
            self.h5f.close()
            print 'closed file'


        

class GuiModel(QtCore.QObject):
    """data_handler -- TraubData object accessing the files.

    ca_cells -- list of cells whose [Ca2+] data is available.

    vm_cells -- list of cells whose Vm data is available.

    spike_cells -- list of cells whose spiketrains are available.
    """
    def __init__(self, *args):
        QtCore.QObject.__init__(self, *args)
        self.data_handler = None
        self.selected_ca = QtGui.QStringListModel()
        self.selected_vm = QtGui.QStringListModel()
        self.selected_spikes = QtGui.QStringListModel()
        self.available_ca = QtGui.QStringListModel()
        self.available_vm = QtGui.QStringListModel()
        self.available_spikes = QtGui.QStringListModel()
        
    def openFile(self, filename):
        """Open a new data file given by filename."""
        self.data_handler = TraubData(filename)
        self.available_spikes.setStringList(self.data_handler.spiketrain_names)
        self.available_vm.setStringList(self.data_handler.Vm_names)
        self.available_ca.setStringList(self.data_handler.Ca_names)
        self.emit(QtCore.SIGNAL('file_loaded(string)'), filename)
        

    # def addTables(self, table_list, table_type):
    #     target = None
    #     source = None
    #     if table_type.lower() == 'ca':
    #         target = self.selected_ca
    #         source = self.available_ca
    #     elif table_type.lower() == 'vm':
    #         target = self.selected_vm
    #         source = self.available_vm
    #     elif table_type.lower() == 'spike':
    #         target = self.selected_spikes
    #         source = self.available_spikes
    #     else:
    #         raise Exception('No such table type: %s' % table_type)
        
    #     for table in table_list:
    #         source.remove(table)
    #         target.append(table)
    #     self.emit(QtCore.SIGNAL('tablelist_updated()'))
        

    # def removeTables(self, table_list, table_type):
    #     target = None
    #     source = None
    #     if table_type.lower() == 'ca':
    #         source = self.selected_ca
    #         target = self.available_ca
    #     elif table_type.lower() == 'vm':
    #         source = self.selected_vm
    #         target = self.available_vm
    #     elif table_type.lower() == 'spike':
    #         source = self.selected_spikes
    #         target = self.available_spikes
    #     else:
    #         raise Exception('No such table type: %s' % table_type)
        
    #     for table in table_list:
    #         source.remove(table)
    #         target.append(table)
    #         target.sort()
    #     self.emit(QtCore.SIGNAL('tablelist_updated()'))
        

class DataVizGui(QtGui.QMainWindow):
    def __init__(self, *args):
        QtGui.QMainWindow.__init__(self, *args)
        self.setDockOptions(self.AllowNestedDocks | self.AllowTabbedDocks | self.ForceTabbedDocks | self.AnimatedDocks)
        self.setDockNestingEnabled(True)
        self.gui_model = GuiModel()
        # self.make_table_lists()
        self.plot = Qwt.QwtPlot()
        self.setCentralWidget(self.plot)
        self.make_actions()
        self.make_menu()

        
    def make_actions(self):
        self.fopen_action = QtGui.QAction(self.tr('&Open'), self)
        self.connect(self.fopen_action, QtCore.SIGNAL('triggered()'), self.popup_fopen)
        self.quit_action = QtGui.QAction(self.tr('&Quit'), self)
        self.connect(self.quit_action, QtCore.SIGNAL('triggered()'), self.do_quit)
        self.show_available_tables_action = QtGui.QAction(self.tr('Show &Available Tables'), self)
        self.show_selected_tables_action = QtGui.QAction(self.tr('Show &Selected Tables'), self)             
        self.show_available_ca_action = QtGui.QAction(self.tr('Show available Ca tables'), self)
        self.show_available_vm_action = QtGui.QAction('Show available Vm tables', self)
        self.show_available_spike_action = QtGui.QAction('Show available spike tables', self)
        self.show_selected_ca_action = QtGui.QAction('Show selected Ca tables', self)
        self.show_selected_vm_action = QtGui.QAction('Show selected Vm tables', self)
        self.show_selected_spike_action = QtGui.QAction('Show selected spike tables', self)
        
        
    def make_menu(self):
        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction(self.fopen_action)
        self.file_menu.addAction(self.quit_action)

        self.view_menu = QtGui.QMenu('&View', self)
        self.view_menu.addAction(self.show_available_tables_action)
        self.view_menu.addAction(self.show_selected_tables_action)
        self.view_menu.addAction(self.show_available_spike_action)
        self.view_menu.addAction(self.show_available_vm_action)
        self.view_menu.addAction(self.show_available_ca_action)
        self.view_menu.addAction(self.show_selected_spike_action)
        self.view_menu.addAction(self.show_selected_vm_action)
        self.view_menu.addAction(self.show_selected_ca_action)
        self.menuBar().addMenu(self.file_menu)
        self.menuBar().addMenu(self.view_menu)
        
    def make_table_lists(self):
        self.selected_ca_view = QtGui.QListView(self)
        self.selected_vm_view = QtGui.QListView(self)
        self.selected_spikes_view = QtGui.QListView(self)
        self.available_ca_view = QtGui.QListView(self)
        self.available_vm_view = QtGui.QListView(self)
        self.available_spikes_view = QtGui.QListView(self)

        self.selected_ca_view.setObjectName('selected_ca')
        self.selected_vm_view.setObjectName('selected_vm')
        self.selected_spikes_view.setObjectName('selected_spikes')
        self.available_ca_view.setObjectName('available_ca')
        self.available_vm_view.setObjectName('available_vm')
        self.available_spikes_view.setObjectName('available_spikes')

        self.selected_ca_view.setModel(self.gui_model.selected_ca)
        self.selected_vm_view.setModel(self.gui_model.selected_vm)
        self.selected_spikes_view.setModel(self.gui_model.selected_spikes)
        self.available_ca_view.setModel(self.gui_model.available_ca)
        self.available_vm_view.setModel(self.gui_model.available_vm)
        self.available_spikes_view.setModel(self.gui_model.available_spikes)
        self.tables_docks = []
        self.list_views = [self.selected_ca_view, self.selected_vm_view,
                 self.selected_spikes_view, self.available_ca_view,
                 self.available_vm_view, self.available_spikes_view]
        for view in self.list_views:
            view.setAcceptDrops(True)
            view.setDragEnabled(True)
            view.setDropIndicatorShown(True)
            view.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
            self.tables_docks.append(QtGui.QDockWidget(view.objectName(), self))
            self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.tables_docks[-1])
        
    def popup_fopen(self):
        file_dialog = QtGui.QFileDialog(self)
        file_dialog.setFileMode(QtGui.QFileDialog.ExistingFile)
        file_dialog.setFilter(self.tr('hdf5(*.h5 *.hd5 *.hdf5);;All files(*)'))
        if file_dialog.exec_():
            filenames = file_dialog.selectedFiles()
            self.gui_model.openFile(str(filenames[0]))

    def do_quit(self):
        print 'Quitting'
        if self.gui_model.data_handler:
            del self.gui_model.data_handler
            print 'Deleted gui model object'
        QtGui.qApp.quit()

            
class __TraubData:
    """Class for visualizing data saved in custom HDF5 files in Traub model simulations.

    spike_train_dict -- a dict which will contain the spike trains in the file as {cell_name: spike_time_array} after loading a data file.

    vm_dict -- a dict which will contain the membrain potential Vm in the file as {cell_name: array_of_Vm_values_over_simulation_time} after loading a file.

    ca_dict -- a dict which will contain the [Ca2+] in the file as {cell_name: array_of_[Ca2+]_values_over_simulation_time} after loading a file.

    plotdt -- plotting time step specified in the data file

    simdt -- simulation time step as specified in the data file.

    timepoints -- an array with simulation time value for each entry in the Vm/Ca arrays.
    
    """
    def __init__(self, filename):
        self.datafilename = filename
        self.spiketrain_dict = {}
        self.vm_dict = {}
        self.ca_dict = {}
        self.spike_matrix = None
        self.simtime = None
        self.simdt = None
        self.plotdt = None
        self.timestamp = None
        self.frame_count = 0
        self.timepoints = []
        self.cell_index_map = {}
        self.index_cell_map = {}
        self._read_data()
        self.draw()
        
    def _read_data(self):
        """Read spiketime serieses, Vm serieses and [Ca2+] serieses from data file."""
        with tables.openFile(self.datafilename) as h5file:
            try:
                self.simtime = h5file.root._v_attrs.simtime
                self.simdt = h5file.root._v_attrs.simdt
                self.plotdt = h5file.root._v_attrs.plotdt
                self.timestamp = h5file.root._v_attrs.timestamp
            except AttributeError, e:
                print e
            try:
                count = 0
                for train in h5file.root.spiketimes:
                    self.spiketrain_dict[train.name] = train[:]
                    cell_name = train.name[:train.name.rfind('_')]
                    self.cell_index_map[cell_name] = count
                    self.index_cell_map[count] = cell_name
            except tables.NoSuchNodeError:
                print 'No node called /spiketimes'
            try:
                for vm_array in h5file.root.Vm:
                    cell_name = vm_array.name[vm_array.name.rfind('_')]
                    self.vm_dict[cell_name] = vm_array[:]
            except tables.NoSuchNodeError:
                print 'No node called /Vm'
            
            try:
                for ca_array in h5file.root.Ca:
                    cell_name = ca_array.name[:ca_array.name.rfind('_')]
                    self.ca_dict[cell_name] = ca_array[:]
            except tables.NoSuchNodeError:
                print 'No node called /Ca'
        print 'SIMULATION DONE ON', self.timestamp, '-- SIMTIME:', self.simtime, ', SIMDT:', self.simdt, ', PLOTDT:', self.plotdt
        
        if len(self.vm_dict) > 0:
            for cell_name, vm_array in self.vm_dict.items():
                self.frame_count = len(vm_array)
                break
        elif self.simtime and self.plotdt:
            self.frame_count = int(self.simtime / self.plotdt + 0.5)
        if not self.simtime:
            self.simtime = 1.0
        self.timepoints = np.linspace(0, self.simtime, self.frame_count)
        if not self.spiketrain_dict:
            pass
            # spike_mat = ll_mat(len(self.timepoints), len(self.spiketrain_dict.keys))
            # for key, value in self.spiketrain_dict.items():
                
        self.vm_axes.set_xlim(self.simtime)
        self.ca_axes.set_xlim(self.simtime)

    def visualize(self):
        """Visualize the data.

        There are three groups -- a set of [Ca2+] recordings, set of
        Vm recordings for the same cells, and spike times for all
        cells.

        I'll display the Vm and [Ca2+] in side by side plots. At the
        same time display all the cells - organized in a grid grouped
        on the basis of depth and celltype.

        """
        # Raster plot for the spike trains: we are skipping the animation for the time being.
        count = 0
        for key, value in self.spiketrain_dict.items():
            self.spike_axes.scatter(self.timepoints, value)
        count = 0
        for key, value in self.vm_dict.items():
            self.vm_axes.plot(self.timepoints, value+count*0.2, label=key)
            count += 1
        count = 0
        for key, value in self.ca_dict.items():
            self.ca_axes.plot(self.timepoints, value+count*0.2, label=key)
            count +=1
        self.draw()


if __name__ == '__main__':
    app = QtGui.QApplication([])
    mainwin = DataVizGui()
    app.lastWindowClosed.connect(mainwin.do_quit)
    mainwin.show()
    app.exec_()
# 
# dataviz.py ends here
