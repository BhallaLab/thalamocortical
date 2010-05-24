# plot.py --- 
# 
# Filename: plot.py
# Description: 
# Author: 
# Maintainer: 
# Created: Wed May 19 08:45:03 2010 (+0530)
# Version: 
# Last-Updated: Mon May 24 11:13:41 2010 (+0530)
#           By: subha
#     Update #: 46
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

import pylab
import config

def plot_data(infilename, outfilename, title=None, xlabel=None, ylabel=None):
    if xlabel is None:
        xlabel = 'time (s)'
    if ylabel is None:
        ylabel = 'Membrane potential (V)'
    if title is None:
        title = infilename
    data = pylab.loadtxt(infilename)
    time_list = pylab.linspace(0,  len(data) * config.plotdt, len(data))
    pylab.plot(time_list, data)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    pylab.title(title)
    pylab.savefig(outfilename)


directory = 'data/20100517'
filenames = [
    'DeepAxoaxonic_0__59.plot',
    'NontuftedRS_0__48.plot',
    'SupAxoaxonic_0__59.plot',
    'SupPyrFRB_0__72.plot',
    'TuftedIB_0__60.plot',
    'DeepBasket_0__59.plot',
    'nRT_0__59.plot',
    'SupBasket_0__59.plot',
    'SupPyrRS_0__72.plot',
    'TuftedRS_0__60.plot',
    'DeepLTS_0__59.plot',
    'SpinyStellate_0__57.plot',
    'SupLTS_0__59.plot',
    'TCR_0__135.plot']
if __name__ == '__main__':
    pylab.hold(False)
    for filename in filenames:
        filename_parts = filename.partition('_')
        outfilename = filename_parts[0] + '.png'
        title = filename_parts[0]
        plot_data(directory + '/' + filename, outfilename, title)
    print 'This is a test change'
# 
# plot.py ends here
