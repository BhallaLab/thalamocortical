#!/usr/bin/env python
# plot.py --- 
# 
# Filename: plot.py
# Description: 
# Author: 
# Maintainer: 
# Created: Wed May 19 08:45:03 2010 (+0530)
# Version: 
# Last-Updated: Wed May 26 03:23:42 2010 (+0530)
#           By: subha
#     Update #: 248
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

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import pylab
import sys
import os
import Tkinter

import config

def destroy(e):
    sys.exit()

class PlotApp:
    """Class to create a little plotter using Tk and matplotlib"""
    def __init__(self, master, directory=None):
        self.current_index = 0
        self.files = []
        if directory is None:
            directory = '.'
        self.set_dir(directory)
        self.jump_step = 10
        self.frame = Tkinter.Frame(master)
        self.frame.pack()
        self.next_button = Tkinter.Button(self.frame, text='NEXT', command=self.plot_next)
        self.prev_button = Tkinter.Button(self.frame, text='PREV', command=self.plot_prev)
        self.jump_label = Tkinter.Label(self.frame, text='Jump steps with PGUP/PGDN')
        self.jump_entry = Tkinter.Entry(self.frame)
        self.jump_entry.insert(0, '10')
        
        self.prev_button.pack(side=Tkinter.LEFT)
        self.next_button.pack(side=Tkinter.LEFT)
        self.jump_entry.pack(side=Tkinter.RIGHT)
        self.jump_label.pack(side=Tkinter.RIGHT)
        self.figure = Figure(figsize=(5,4), dpi=100)
        self.axes = self.figure.add_subplot(111)
        # self.subplot.hold(False)
        self.canvas = FigureCanvasTkAgg(self.figure, master)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        self.toolbar = NavigationToolbar2TkAgg( self.canvas, master )
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)


        self.xlabel = 'Time (s)'
        self.ylabel = 'Membrane potential (V)'
        self.title = None
        self.do_plot(self.files[self.current_index])

    def set_dir(self, path):
        self.files = [path + '/' + filename for  filename in os.listdir(path)]
        self.files.sort()
        self.current_index = 0

    def plot_next(self, step=1):
        if not self.files:
            print 'No files to plot'
            return
        
        self.current_index = (self.current_index + step) % len(self.files)
        self.do_plot(self.files[self.current_index])

    def plot_prev(self, step=1):
        if not self.files:
            print 'No files to plot'
            return
        
        self.current_index = (self.current_index - step) % len(self.files)
        self.do_plot(self.files[self.current_index])

    def do_plot(self, filename):
        # print 'Plotting', filename
        try:
            data = pylab.loadtxt(filename)
            time_list = pylab.linspace(0, len(data) * config.plotdt, len(data))
            self.axes.cla()
            self.title = self.files[self.current_index]
            self.axes.set_xlabel(self.xlabel)
            self.axes.set_ylabel(self.ylabel)
            self.axes.set_title(self.title)
            self.axes.plot(time_list, data)
            self.canvas.draw()
            # print 'Title:', self.axes.get_title()
        except Exception, e:
            print e
            # print 'FILES:', self.files
    
    def handle_leftright_keys(self, event):
        if event.keysym == 'Left':
            self.plot_prev()
        elif event.keysym == 'Right':
            self.plot_next()
        if event.keysym == 'Prior':
            self.plot_prev(self.jump_step)
        elif event.keysym == 'Next':
            self.plot_next(self.jump_step)
            

if __name__ == '__main__':
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        directory = '.'
    root = Tkinter.Tk()
    root.wm_title('Cortical model plots using matplotlib and Tk')
    app = PlotApp(root, directory)
    root.bind('<Destroy>', destroy)
    root.bind('<KeyPress>', app.handle_leftright_keys)
    root.mainloop()


# 
# plot.py ends here
