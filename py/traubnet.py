# traubnet.py --- 
# 
# Filename: traubnet.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Aug 10 15:45:05 2010 (+0530)
# Version: 
# Last-Updated: Wed Sep 15 08:29:45 2010 (+0530)
#           By: subha
#     Update #: 760
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# This replaces TraubNet class in connection.py. I am switching from
# strings and maps to networkx graphs and attributes of the nodes and
# the edges.
# 
# 

# Change log:
# 
# 2010-08-10 15:48:26 (+0530) - initiated development.
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

from datetime import datetime
import sys
import os
import csv
import numpy
import matplotlib
import allowedcomp

# matplotlib.use('SVG')

import matplotlib.pyplot as plt
has_mayavi = False
# try:
#     from enthought.mayavi import mlab
#     has_mayavi = True
#    mlab.options.offscreen = True #-- this causes mayavi to crash
# known bug: https://svn.enthought.com/enthought/ticket/1824
# except ImportError:
#     has_mayavi = False
#     print 'Mayavi Python Module not available on this system. No 3-D visualization'
    

import networkx as nx

import config
import synapse

class TraubNet(object):
    """
    Store network information in graphs.
    
    """
    def __init__(self, 
                 celltypes_file='nx_celltype_graph.edgelist', 
                 cells_file='nx_cell_graph.edgelist', format='edgelist', scale=1.0):
        """
        celltypes_file -- file containing celltype-celltype connectivity graph
        
        cells_file -- file containing cell-cell connectivity graph

        format -- string representation of file format of the celltypes_file and cells_file.

        """
        self.__celltype_graph = self._read_celltype_graph(celltypes_file, format=format)
        if not self.__celltype_graph:
            self.__celltype_graph = self._make_celltype_graph('connmatrix.txt', 'cells.txt', scale=scale)
        self.__cell_graph = self._read_cell_graph(cells_file, format=format)
        if not self.__cell_graph:
            self.__cell_graph = self._make_cell_graph()
        start = datetime.now()
        print nx.info(self.__cell_graph)
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Computed Graph info in: %g' % (delta.seconds + 1e-6 * delta.microseconds))

    def _make_celltype_graph(self, connmatrix_file, cellcount_file, scale=1.0):        
        """
        connmatrix_file -- csv file containig the connection matrix.
        
        The first row in the csv file should be the column headers - which
        are the cell types. The connection matrix itself is a square
        matrix with entry[i][j] specifying the number of presynaptic
        cells per postsynaptic cell, where the presynaptic cells are
        of type header[i] and the postsynaptic cell is of type
        header[j]

        cellcount_file -- csv file containing the number of instances
        of each celltype.

        scale -- reduce/increase the number of cells of each type by this multiple

        """
        start = datetime.now()
        celltype_graph = None
        cellcount_dict = {}
        with open(cellcount_file, 'r') as popfile:
            reader = csv.reader(popfile)
            for line in reader:
                if len(line) > 0:
                    cellcount_dict[line[0]] = int(line[1])
        celltype_graph = nx.DiGraph()
        index = 0
        for key, value in cellcount_dict.items():
            print key, value
            value = int(numpy.round(value*scale))
            celltype_graph.add_node(key, count=value, index=index)
            index += 1
        reader = csv.reader(file(connmatrix_file))
        header = reader.next()
        row = 0
        for line in reader:
            if len(line) <= 0:
                continue
            pre = header[row]
            row += 1
            col = 0
            for entry in line:
                post = header[col]
                col += 1
                value = int(numpy.round(int(entry)*scale))
                if value == 0:
                    config.LOGGER.debug('No edge between %s and %s' % (pre, post))
                    continue
                celltype_graph.add_edge(pre, post, weight=value)

                try:
                    allowed_comps = allowedcomp.ALLOWED_COMP[pre][post]
                except KeyError:
                    allowed_comps = []
                celltype_graph[pre][post]['ps_comps'] = str(allowed_comps)

                if pre != 'nRT':
                    try:
                        tau_gaba_fast = synapse.TAU_GABA[pre][post]                            
                        celltype_graph[pre][post]['tau_gaba_fast'] = tau_gaba_fast
                    except KeyError:
                        config.LOGGER.info('No tau_gaba for synapse between %s and %s' % (pre, post))
                else:
                    try:
                        tau_gaba_fast = synapse.TAU_GABA_FAST[pre][post]
                        celltype_graph[pre][post]['tau_gaba_fast'] = tau_gaba_fast
                    except KeyError:
                        config.LOGGER.info('No tau_gaba for synapse between %s and %s' % (pre, post))                            
                    try:
                        tau_gaba_slow = synapse.TAU_GABA_SLOW[pre][post]
                        celltype_graph[pre][post]['tau_gaba_slow'] = tau_gaba_slow
                    except KeyError:
                        config.LOGGER.info('No tau_gaba for synapse between %s and %s' % (pre, post))
                try:
                    tau_ampa = synapse.TAU_AMPA[pre][post]
                    celltype_graph[pre][post]['tau_ampa'] = tau_ampa
                except KeyError:
                    config.LOGGER.info('No tau_ampa for synapse between %s and %s' % (pre, post))
                try:
                    tau_nmda = synapse.TAU_NMDA[pre][post]
                    celltype_graph[pre][post]['tau_nmda'] = tau_nmda
                except KeyError:
                    config.LOGGER.info('No tau_gaba for synapse between %s and %s' % (pre, post))
                try:
                    gbar_gaba = synapse.G_GABA[pre][post]
                    celltype_graph[pre][post]['gbar_gaba'] = gbar_gaba
                except KeyError:
                    config.LOGGER.info('No gbar_gaba for synapse between %s and %s' % (pre, post))
                try:
                    gbar_ampa = synapse.G_AMPA[pre][post]
                    celltype_graph[pre][post]['gbar_ampa'] = gbar_ampa
                except KeyError:
                    config.LOGGER.info('No gbar_ampa for synapse between %s and %s' % (pre, post))
                try:
                    gbar_nmda = synapse.G_NMDA[pre][post]
                    celltype_graph[pre][post]['gbar_nmda'] = gbar_nmda
                except KeyError:
                    config.LOGGER.info('No gbar_nmda for synapse between %s and %s' % (pre, post))
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Generated celltype_graph in %g s' % (delta.seconds + delta.microseconds * 1e-6))
        return celltype_graph
        
    def _read_celltype_graph(self, 
                             celltypes_file, 
                             format='gml'):
        """
        Read celltype-celltype connectivity graph from file.

        celltypes_file -- the path of the file containing
        the graph.
        
        format -- format of the file. allowed values: gml, graphml, edgelist, pickle, yaml.

        """
        start = datetime.now()
        celltype_graph = None
        try:
            if format == 'gml':
                celltype_graph = nx.read_gml(celltypes_file)
            elif format == 'edgelist':
                celltype_graph = nx.read_edgelist(celltypes_file)
            elif format == 'graphml':
                celltype_graph = nx.read_graphml(celltypes_file)
            elif format == 'pickle':
                celltype_graph = nx.read_gpickle(celltypes_file)
            elif format == 'yaml':
                celltype_graph = nx.read_yaml(celltypes_file)
            else:
                print 'Unrecognized format %s' % (format)
        except Exception, e:
            print e
        if celltype_graph is not None:
            if not ('doc' in celltype_graph):
                celltype_graph.graph['doc'] = 'Celltype-based connectivity data. \
count of node *n* is the number of cells of type *n* \
that are present in the model. weight of edge (a, b) \
is the number of cells of type *a* that connect to \
each cell of type *b*.'
            if not celltype_graph.name:
                celltype_graph.name = 'CellTypeGraph'
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Read celltype_graph from file %s of format %s in %g s' 
                                     % (celltypes_file, format, delta.seconds + 1e-6 * delta.microseconds))
        return celltype_graph

    def plot_celltype_graph(self, filename='celltype_graph.png'):
        """Display the celltype connectivity graph 

        """
        try:
            pos = nx.graphviz_layout(self.__celltype_graph)
        except Exception, e:
            print e
            pos = nx.spring_layout(self.__celltype_graph)
        node_size_list = [self.__celltype_graph.node[vertex]['count'] * 10 for vertex in self.__celltype_graph]
        edge_weights = [edata['weight'] for u, v, edata in self.__celltype_graph.edges(data=True)]
        plt.figure(1)
        nx.draw(self.__celltype_graph,
                pos,
                alpha=0.4,
                with_labels=True,
                node_size = node_size_list,
                edge_color=edge_weights,
                edge_cmap=plt.cm.jet,
                edge_vmin=min(edge_weights),
                edge_vmax=max(edge_weights))
        #plt.show()
        plt.savefig(filename)
        print 'Saved figure in', filename

    def plot_celltype_graph_3d(self, filename='celltypes_graph_3d.png'):
        """Some eyecandi useful for presentations."""
        if not has_mayavi:
            return
        numeric_graph = nx.convert_node_labels_to_integers(self.__celltype_graph)
        pos = nx.spring_layout(numeric_graph, dim=3)        
        xyz = numpy.array([pos[v] for v in numeric_graph])
        scalars = [self.__celltype_graph.node[vertex]['count'] for vertex in self.__celltype_graph]
        fig = mlab.figure(1, bgcolor=(0, 0, 0))
        mlab.clf()
        points = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                               scalars,
                               scale_factor=0.1,
                               scale_mode='none',
                               colormap='summer',
                               opacity=0.4,
                               transparent=True,
                               resolution=20)
        points.mlab_source.dataset.lines = numpy.array(numeric_graph.edges())
        points.mlab_source.update()
        # mlab.pipeline.surface(points, color=(1,1,1),
        #                       representation='wireframe',
        #                       line_width=2,
        #                       name='synapses')
        tube = mlab.pipeline.tube(points, tube_radius=0.01)
        mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
        mlab.savefig(filename, size=(1280, 800), figure=fig)
        print 'Mayavi celltype graph saved in', filename
        mlab.show()

    def save_celltype_graph(self, filename='celltype_conn.gml', format='gml'):
        """
        Save the celltype-to-celltype connectivity information in a file.
        
        filename -- path of the file to be saved.

        format -- format to save in. Using GML as GraphML support is
        not complete in NetworkX.  

        """
        start = datetime.now()
        if format == 'gml':
            nx.write_gml(self.__celltype_graph, filename)
        elif format == 'yaml':
            nx.write_yaml(self.__celltype_graph, filename)
        elif format == 'graphml':
            nx.write_graphml(self.__celltype_graph, filename)
        elif format == 'edgelist':
            nx.write_edgelist(self.__celltype_graph, filename)
        elif format == 'pickle':
            nx.write_gpickle(self.__celltype_graph, filename)
        else:
            raise Exception('Supported formats: gml, graphml, yaml. Received: %s' % (format))
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Saved celltype_graph in file %s of format %s in %g s' 
                                     % (filename, format, delta.seconds + delta.microseconds * 1e-6))
        print 'Saved celltype connectivity graph in', filename

    def _read_cell_graph(self, filename, format):
        """Load the cell-to-cell connectivity graph from a
        file. 

        Returns None if any error happens.
        """
        cell_graph = None        
        if filename:
            try:
                start = datetime.now()
                if format == 'gml':
                    cell_graph = nx.read_gml(filename)
                elif format == 'pickle':
                    cell_graph = nx.read_gpickle(filename)
                elif format == 'edgelist':
                    cell_graph = nx.read_edgelist(filename)
                elif format == 'yaml':
                    cell_graph = nx.read_yaml(filename)
                elif format == 'graphml':
                    cell_graph = cell_graph = nx.read_graphml(filename)
                else:
                    print 'Unrecognized format:', format
                end = datetime.now()
                delta = end - start
                config.BENCHMARK_LOGGER.info('Read cell_graph from file %s of format %s in %g s' 
                                             % (filename, format, delta.seconds + 1e-6 * delta.microseconds))
            except Exception, e:
                print e
        return cell_graph

    def _make_cell_graph(self):
        """Expand the celltype-to-celltype connectivity information
        and make a graph representing the network of the individual
        cells.

        Each cell is identified by the string {celltype}_{index}

        """
        print 'Creating the cell_graph from scratch'
        start = datetime.now()
        cell_graph = nx.MultiDiGraph()
        cell_graph.name = 'CellGraph'
        for celltype in self.__celltype_graph:
            for index in range(self.__celltype_graph.node[celltype]['count']):
                cell_graph.add_node('%s_%d' % (celltype, index), type_index=self.__celltype_graph.node[celltype]['index'])

        for pre, post, edata in self.__celltype_graph.edges(data=True):
            pre_count = self.__celltype_graph.node[pre]['count']
            post_count = self.__celltype_graph.node[post]['count']            
            pre_post_ratio = edata['weight']
            ps_comps = numpy.array(eval(edata['ps_comps']), dtype=int)
            if (pre_post_ratio == 0) or (len(ps_comps) == 0):
                continue
            # randint returns unifrom random integers in [low, high)
            # interval. i-th row of pre_indices = list of indices of
            # presynaptic cells of type pre connected to i-th cell of
            # type post.
            pre_indices = numpy.random.randint(low=0, high=pre_count, size=(post_count, pre_post_ratio)) 
            # List of indices in ps_comps for all post-synaptic cells
            # for all presynaptic cells.
            # The i-th row of post_comp_indices will have the list of
            # indices of the post synaptic compartments for each
            # presynaptic cell of i-th post cell. Thus, if
            # ps_comp_indices[i][j] = k, then the j-th pre cell for i-th
            # post cell will create a synapse on k-th entry in
            # the ps_comps list for the pre-post edge.
            post_comp_indices = numpy.random.randint(low=0, high=len(ps_comps), size=(post_count, pre_post_ratio))
            
            for ii in range(post_count):
                post_cell_name = '%s_%d' % (post, ii)
                post_comps = ps_comps[post_comp_indices[ii]]
                pre_indices_for_post = pre_indices[ii]
                for jj in range(pre_post_ratio):
                    pre_cell_name = '%s_%d' % (pre, pre_indices[ii][jj])
                    cell_graph.add_edge(pre_cell_name, post_cell_name, ps_comp=post_comps[jj])
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Built cell_graph programmatically - time: %g s' % (delta.seconds + 1e-6 * delta.microseconds))
        return cell_graph

    def plot_cell_graph(self):
        """Display the cell-to-cell connection graph.

        """
        plt.figure(3)
        nx.draw(self.__cell_graph,
                alpha=0.4,
                with_labels=False,
                node_color=[self.__cell_graph.node[vertex]['type_index'] for vertex in self.__cell_graph],
                cmap=plt.cm.jet,
                vmin=0,
                vmax=len(self.__celltype_graph),
                edge_color='b')
        #plt.show()
        plt.savefig('cell_graph')
        print 'Saved 2D cell graph in', 'cell_graph'

    def plot_cell_graph_3d(self, filename=None):
        """Some eyecandi useful for presentations."""
        if not has_mayavi:
            return
        numeric_graph = nx.convert_node_labels_to_integers(self.__cell_graph)
        pos = nx.spring_layout(numeric_graph, dim=3)        
        xyz = numpy.array([pos[v] for v in numeric_graph])
        celltypes = self.__celltype_graph.nodes()

        scalars = []
        for cell in self.__cell_graph:
            celltype = cell.split('_')[0]
            scalars.append(celltypes.index(celltype))
        scalars = numpy.array(scalars)
        fig = mlab.figure(4, size=(1280, 800), bgcolor=(0, 0, 0))
        mlab.clf()
        # print xyz
        # print xyz[:,0], xyz[:,1], xyz[:,2]
        points = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                               scalars,
                               scale_factor=0.01,
                               scale_mode='none',
                               colormap='autumn',
                               opacity=0.5,
                               resolution=20,
                               vmin = 0.0,
                               vmax=len(celltypes),
                               transparent=True,
                               )
        points.mlab_source.dataset.lines = numpy.array(numeric_graph.edges())
        tube = mlab.pipeline.tube(points, tube_radius=0.001)
        mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
        camera = mlab.view()
        mlab.view(camera[0], camera[1], camera[2] - 1.0, camera[3])
        print 'Azimuth, Elevation, Distance, Focal point', camera
	# mlab.move(3) # Move the camera forward
        if filename is None:
            filename = 'cell_graph_3d.png'
        mlab.savefig(filename, size=(1280,800), figure=fig)
        print 'Mayavi image saved in', filename
        #mlab.show()
        

    def save_cell_graph(self, filename='cell_graph.gml', format='gml'):
        """Save the cell to cell connectivity graph in a file.

        """
        start = datetime.now()
        if format == 'gml':
            nx.write_gml(self.__cell_graph, filename)
        elif format == 'pickle':
            nx.write_gpickle(self.__cell_graph, filename)
        elif format == 'yaml':
            nx.write_yaml(self.__cell_graph, filename)
        elif format == 'graphml':
            nx.write_graphml(self.__cell_graph, filename)
        elif format == 'edgelist':
            nx.write_edgelist(self.__cell_graph, filename)
        else:
            raise Exception('Not supported: %s' % (format))
        end = datetime.now()
        delta = end - start
        config.BENCHMARK_LOGGER.info('Saved cell graph in file %s of type %s in %g s' %(filename, format, delta.seconds + 1e-6 * delta.microseconds))
        print 'Saved cell-to-cell connectivity data in', filename
    
def test(args=None):
    """the first argument specifies the graph file format. Default is GML

    second argument is the scale factor. All the cell and edge counts
    in the celltype graph are multiplied by this factor.

    third argument was used earlier to save the cell-graph image in a
    file.
    
    """
    if len(args) > 1:
        format = args[1]
    else:
        format = 'gml'
    if len(args) > 2:
        scale = float(args[2])
    else:
        scale = 1.0
    if len(args) > 3:
        filename = args[3]
    else:
        filename = 'cell_graph.png'
    celltype_graph_file = 'nx_celltype_graph.' + format
    cell_graph_file = 'nx_cell_graph.' + format
    net = TraubNet(celltype_graph_file, cell_graph_file, format=format, scale=scale)    
    # net.plot_celltype_graph()
    # net.plot_celltype_graph_3d()
    # net.save_celltype_graph(filename=celltype_graph_file, format=format)
    # net.plot_cell_graph()
    # net.plot_cell_graph_3d()
    # net.save_cell_graph(cell_graph_file, format=format)

if __name__ == '__main__':
    """the first argument specifies the graph file format. Default is GML

    second argument is the scale factor. All the cell and edge counts
    in the celltype graph are multiplied by this factor.

    third argument was used earlier to save the cell-graph image in a
    file.
    
    """
    print sys.argv
    test(sys.argv)

# 
# traubnet.py ends here
