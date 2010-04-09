# test_copy.py --- 
# 
# Filename: test_copy.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Apr  8 19:01:04 2010 (+0530)
# Version: 
# Last-Updated: Fri Apr  9 16:17:11 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 19
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

import config
import moose
from spinystellate import SpinyStellate

def test_copy():
    cellcount = 100
    proto = SpinyStellate.prototype
    cells = []
    config.BENCHMARK_LOGGER.info('TEST1: Starting making cell copies.')
    for ii in range(cellcount):
        cell = moose.Cell(proto, 'cell' + str(ii))
        cells.append(cell)
    config.BENCHMARK_LOGGER.info('TEST1: Finished making cell copies.')
    config.BENCHMARK_LOGGER.info('TEST2.A: Starting making TraubCell copies.')
    ss = []
    for ii in range(cellcount):
        cell = SpinyStellate(proto, 'ss' + str(ii))
        ss.append(cell)
    config.BENCHMARK_LOGGER.info('TEST2.A: Finished making TraubCell copies.')
    config.BENCHMARK_LOGGER.info('TEST2.B: Starting making TraubCells fresh.')
    ss = []
    for ii in range(cellcount):
        cell = SpinyStellate('freshss' + str(ii))
        ss.append(cell)
    config.BENCHMARK_LOGGER.info('TEST2.B: Finished making TraubCells fresh.')        
if __name__ == '__main__':
    test_copy()

# 
# test_copy.py ends here
