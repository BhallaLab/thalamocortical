#!/usr/bin/env python

# Filename: dataviz.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Wed Dec 15 10:16:41 2010 (+0530)
# Version: 
# Last-Updated: Mon Feb 14 20:22:10 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 700
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
        

class DataVizGui(QtGui.QMainWindow):
    def __init__(self, *args):
        QtGui.QMainWindow.__init__(self, *args)
        self.setDockOptions(self.AllowNestedDocks | self.AllowTabbedDocks | self.ForceTabbedDocks | self.AnimatedDocks)
        self.setDockNestingEnabled(True)
        self.gui_model = GuiModel()
        self.make_table_lists()
        self.plot_panel = QtGui.QTabWidget(self)
        self.setCentralWidget(self.plot_panel)
        self.plot = Qwt.QwtPlot(self.plot_panel)
        self.plot_panel.addTab(self.plot, 'plot')
        self.make_actions()
        self.make_menu()
        self.toolbar = self.make_toolbar()
        
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
        self.select_by_regex_action = QtGui.QAction('Select by regex', self)
        
        
        
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
            self.tables_docks[-1].setWidget(view)
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
        QtGui.qApp.quit()

            


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    QtGui.qApp = app
    mainwin = DataVizGui()
    app.lastWindowClosed.connect(mainwin.do_quit)
    mainwin.show()
    app.exec_()
# 
# dataviz.py ends here
