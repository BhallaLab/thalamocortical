# test_mayavi.py --- 
# 
# Filename: test_mayavi.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Aug 28 15:54:27 2010 (+0530)
# Version: 
# Last-Updated: Sat Aug 28 16:42:41 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 60
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

import networkx as nx
import numpy as np
# from enthought.etsconfig.api import ETSConfig
# ETSConfig.toolkit = 'null'
from enthought.mayavi import mlab
 
# mlab.options.offscreen = True

if __name__ == '__main__':
    g = nx.DiGraph()
    g.add_edges_from([(0, 1), (0, 2), (0, 3), 
                      (1, 2), (1,3),
                      (2, 3)])
    pos = nx.spring_layout(g, dim=3)
    xyz = np.array([pos[v] for v in g])
    fig = mlab.figure(1, bgcolor=(0,0,0))
    mlab.clf()
    points = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                           [v for v in g],
                           scale_factor=0.1,
                           scale_mode='none',
                           colormap='summer',
                           resolution=20)
    print points
    edge_start = []
    edge_end = []
    for edge in g.edges():
        edge_start.append(pos[edge[0]])
        edge_end.append(pos[edge[1]])

    edge_start = np.array(edge_start)
    edge_end = np.array(edge_end)
    components = edge_end - edge_start
    mlab.quiver3d(edge_start[:,0], edge_start[:,1], edge_start[:,2], components[:,0], components[:,1], components[:,2],
                  scale_factor=0.1
                  )
    mlab.show()
    print xyz
    


# 
# test_mayavi.py ends here
