# test_viz.py --- 
# 
# Filename: test_viz.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Copyright (C) 2010 Subhasis Ray, all rights reserved.
# Created: Fri Dec 31 17:43:01 2010 (+0530)
# Version: 
# Last-Updated: Fri Dec 31 17:47:48 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 15
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

# Code:

import os
import sys
import time

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from PyQt4 import QtCore, QtGui


ITERS = 50
YMAX = 3500
XMAX = 1000
class BlitQT(FigureCanvas):

    def __init__(self):
        FigureCanvas.__init__(self, Figure())

        self.ax = self.figure.add_subplot(111)
        self.draw()

        self.old_size = self.ax.bbox.width, self.ax.bbox.height
        self.ax_background = self.copy_from_bbox(self.ax.bbox)
        self.cnt = 0
        self.data = np.random.rand(YMAX, XMAX)
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
        self.im.set_data(np.random.rand(YMAX, XMAX))
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
    widget = BlitQT()
    widget.show()

    sys.exit(app.exec_())
    
# 
# dataviz.py ends here



# 
# test_viz.py ends here
