#!/usr/bin/env python

# Filename: dataviz.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Wed Dec 15 10:16:41 2010 (+0530)
# Version: 
# Last-Updated: Tue Dec 21 16:57:48 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 283
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
# simulation code and the data analysis code.
#    Hence forking this away into code for analyzing older data.

# Code:

import os
import sys
import time

import numpy as np
from pysparse.sparse.pysparseMatrix import PysparseMatrix
import tables
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from PyQt4 import QtCore, QtGui


ITERS = 1000

class TraubDataViz(FigureCanvas):
    """Class for visualizing data saved in custom HDF5 files in Traub model simulations."""
    def __init__(self, filename):
        FigureCanvas.__init__(self, Figure())
        self.spike_axes = self.figure.add_subplot(131)
        self.spike_axes.set_title('Spike trains')
        self.vm_axes =  self.figure.add_subplot(132)
        self.vm_axes.set_title('Vm')
        self.ca_axes = self.figure.add_subplot(133)
        self.ca_axes.set_title('[Ca2+]')
        self.spike_axes_bg = self.copy_from_bbox(self.spike_axes.bbox)
        self.vm_axes_bg = self.copy_from_bbox(self.vm_axes.bbox)
        self.ca_axes_bg = self.copy_from_bbox(self.ca_axes.bbox)
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
                    cell_name = train.name[:train.name.rfind('_')]
                    self.spiketrain_dict[cell_name] = train[:]
                    self.cell_index_map[cell_name] = count
                    self.index_cell_map[count] = cell_name
                    count += 1
            except tables.NoSuchNodeError:
                print 'No node called /spiketimes'
            try:
                for vm_array in h5file.root.Vm:
                    cell_name = vm_array.name[:vm_array.name.rfind('_')]
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
        
        if self.vm_dict:
            for cell_name, vm_array in self.vm_dict.items():
                self.frame_count = len(vm_array)
                break
        elif self.simtime and self.plotdt:
            self.frame_count = int(self.simtime / self.plotdt + 0.5)
        if not self.simtime:
            self.simtime = 1.0
        self.timepoints = np.linspace(0, self.simtime, self.frame_count)
        if self.spiketrain_dict:
            spike_mat = PysparseMatrix(nrow=len(self.spiketrain_dict.keys()), ncol=len(self.timepoints))
            for index in range(len(self.index_cell_map.keys())):
                cell_name = self.index_cell_map[index]
                try:
                    spiketrain = self.spiketrain_dict[cell_name]
                    spike_mat.put(1.0,  np.array([index] * len(spiketrain), dtype='int32'), np.cast['int32'](spiketrain/self.plotdt + 0.5))
                except KeyError:
                    print 'No cell corresponding to index', index
            self.spike_matrix = spike_mat.getNumpyArray()
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
        if self.spiketrain_dict:
            for cellname, spiketrain in self.spiketrain_dict.items():
                y_values = np.array([self.cell_index_map[cellname]] * len(spiketrain))
                self.spike_axes.plot(spiketrain, y_values, '.')
                # self.spike_axes.scatter(spiketrain, y_values, s=1, c=y_values, vmin=0.0, vmax=len(y_values))
        print 'finished spike plot'
        count = 0
        # for key, value in self.vm_dict.items():
        #     self.vm_axes.plot(self.timepoints, value+count*0.2, label=key)
        #     count += 1
        # count = 0
        # for key, value in self.ca_dict.items():
        #     self.ca_axes.plot(self.timepoints, value+count*0.2, label=key)
        #     count +=1
        self.draw()
        self.figure.savefig(self.datafilename + '_.png')

# The BlitQt is modified from the matplotlib example -- left here as a
# reference when modifying code in future.
class BlitQT(FigureCanvas):

    def __init__(self):
        FigureCanvas.__init__(self, Figure())

        self.ax = self.figure.add_subplot(111)
        self.draw()

        self.old_size = self.ax.bbox.width, self.ax.bbox.height
        self.ax_background = self.copy_from_bbox(self.ax.bbox)
        self.cnt = 0
        self.data = np.random.rand(350, 100)
        self.im = self.ax.imshow(self.data, aspect='auto', interpolation='nearest', animated=True)
        self.cbar = self.figure.colorbar(self.im)

        self.draw()
        self.tstart = time.time()
        self.startTimer(10)

    def timerEvent(self, evt):
        # The following lines from the example were causing the screen to become blank.
        # Hence I commented them out. -- Subha
        
        # current_size = self.ax.bbox.width, self.ax.bbox.height
        # if self.old_size != current_size:
        #     self.old_size = current_size
        #     self.ax.clear()
        #     # self.ax.grid()
        #     self.draw()
        #     self.ax_background = self.copy_from_bbox(self.ax.bbox)

        self.restore_region(self.ax_background, bbox=self.ax.bbox)

        # update the data
        self.im.set_data(np.random.rand(350, 100))
        # just draw the animated artist -- not needed when updating everything -- Subha
        # self.ax.draw_artist(self.sin_line)
        # self.ax.draw_artist(self.cos_line)
        # just redraw the axes rectangle
        self.blit(self.ax.bbox)
        if self.cnt==ITERS:
            # print the timing info and quit
            print 'FPS:' , ITERS/(time.time()-self.tstart)
            sys.exit()
        else:
            self.cnt += 1
            self.draw()
            self.figure.savefig('frame_%d.png' % (self.cnt))

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    if len(sys.argv) > 1:
        widget = TraubDataViz(sys.argv[1])
        widget.visualize()
        widget.show()

    sys.exit(app.exec_())
    
# 
# dataviz.py ends here
