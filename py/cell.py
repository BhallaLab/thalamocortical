# cell.py --- 
# 
# Filename: cell.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Jul 24 10:04:47 2009 (+0530)
# Version: 
# Last-Updated: Wed Apr 21 10:25:32 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 401
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# This is an extension of Cell class - to add some utility 
# functions for debugging. All cell types should derive from this.
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
from collections import defaultdict

import moose
import config
import pymoose

from nachans import *
from kchans import *
from cachans import *
from archan import *
from capool import *
from compartment import MyCompartment

def init_channel_lib():
    """Initialize the prototype channels in library"""
    if not config.channel_lib:
        config.LOGGER.debug('* Generating channel prototypes in /library')
        for channel_name in config.channel_name_list:
            channel_class = eval(channel_name)
            channel = channel_class(channel_name, config.lib)
            config.channel_lib[channel_name] = channel
            config.LOGGER.debug( '* Created %s' % (channel.path))
        config.channel_lib['SpikeGen'] = moose.SpikeGen('spike', config.lib)
    return config.channel_lib

def nameindex(comp):
    """Utility function to sort by index in the compartment name"""
    if comp is None:
        return -1
    pos = comp.name.rfind('_')
    if pos >= 0:
        index = int(comp.name[pos+1:])
        return index
    else:
        return -1

def get_comp(cell, index):
    """Return a wrapper over compartment specified by index. None if
    no such compartment exists."""
    path = cell.path + '/comp_' + str(index)
    if config.context.exists(path):
        return MyCompartment(path)
    else:
        return None


class TraubCell(moose.Cell):
    channel_lib = init_channel_lib()
                     
    def __init__(self, *args):
        print args
        moose.Cell.__init__(self, *args)
        
    
    # Dynamic access to a compartment by index.  It mimics a python
    # list 'comp' via underlying function call to get_comp(cell,
    # index)
    comp = moose.listproperty(get_comp)

    @property
    def soma(self):
        return get_comp(self, 1)

    def pfile_name(self):
        """Each cell type subclass should implement this"""
        raise NotImplementedError, "function pfile_name not implemented"

    @classmethod
    def read_proto(cls, filename, cellname, params=None):
        """Read a prototype cell from .p file into library.  

        Each cell type class should initialize its prototype with a
        call to this function. with something like this within the
        class declaration:

        prototype = TraubCell.read_proto("MyCellType.p", "MyClassName")

        filename -- path(relative/absolute) of the cell prototype file.

        cellname -- path of the cell to be Created

        params -- if specified, channels in /library are adjusted with
        the parameters specified in this (via a call to
        adjust_chanlib).

        """
        config.LOGGER.debug('Reading proto:%s' % (filename))
        if params is not None:
            TraubCell.adjust_chanlib(params)
        ret = None
        cellpath = config.lib.path + '/' + cellname
        if not config.context.exists(cellpath):
            config.LOGGER.debug(__name__ + ' reading cell: ' + cellpath)
            for handler in config.LOGGER.handlers:
                handler.flush()
            config.context.readCell(filename, cellpath)
        else:
            config.LOGGER.debug(__name__ + ' cell exists: ' + cellpath)
	config.LOGGER.debug('Returning cell %s' % (cellpath))
        for handler in config.LOGGER.handlers:
            handler.flush()
        return moose.Cell(cellpath)

    @classmethod
    def adjust_chanlib(cls, chan_params):
        """Set the properties of prototype channels in /library to fit
        the channel properties of this cell type.

        chan_params -- dict containing the channel parameters. The
        following string keys should be there with float values:

        ENa
        EK
        EAR
        ECa
        TauCa
        
        """
        config.LOGGER.debug('Adjusting channel properties.')
        for key, channel in init_channel_lib().items():
            if isinstance(channel, KChannel):
                channel.Ek = chan_params['EK']
            elif isinstance(channel, NaChannel):
                channel.Ek = chan_params['ENa']
            elif isinstance(channel, CaChannel):
                channel.Ek = chan_params['ECa']
            elif isinstance(channel, AR):
                channel.Ek = chan_params['EAR']
                try:
                    channel.X = chan_params['X_AR']
                except KeyError:
                    channel.X = 0.25
            elif isinstance(channel, CaPool):
                channel.tau = chan_params['TauCa']

    

    def _ca_tau(self):
        raise NotImplementedError("You must set tau for [Ca2+] decay in the method _ca_tau() in subclass.")

    def _setup_passive(self):
        raise NotImplementedError("You must define _setup_passive to set the passive membrane properties and other post-readcell tweakings.")

    def _setup_channels(self):
        raise NotImplementedError("You must define setup_channels to set the channel reversal potential and other post-readcell tweakings.")

    def _topology(self):
        raise NotImplementedError("You must define cell topology in the method _topology() in subclass.")

    def has_cycle(self, comp=None):
        if comp is None:
            comp = self.soma
        comp._visited = True
        ret = False
        for item in comp.raxial_list:
            if hasattr(item, '_visited') and item._visited:
                config.LOGGER.warning('Cycle between: %s and %s.' % (comp.path, item.path))
                return True
            ret = ret or has_cycle(item)
        return ret

    def dump_cell(self, file_path):
        dump_file = open(file_path, "w")
        dump_file.write("comp,len,dia,sarea,xarea,Em,Cm,Rm,Ra")
        dump_file.write(",e_ar,gbar_ar")
        dump_file.write(",tau_cad,beta_cad")
        dump_file.write(",e_cal,gbar_cal")
        dump_file.write(",e_cat,gbar_cat")
        dump_file.write(",e_cat_a,gbar_cat_a")
        dump_file.write(",e_k2,gbar_k2")
        dump_file.write(",e_ka,gbar_ka")
        dump_file.write(",e_ka_ib,gbar_ka_ib")
        dump_file.write(",e_kahp,gbar_kahp")
        dump_file.write(",e_kahp_deeppyr,gbar_kahp_deeppyr")
        dump_file.write(",e_kahp_slower,gbar_kahp_slower")
        dump_file.write(",e_kc,gbar_kc")
        dump_file.write(",e_kc_fast,gbar_kc_fast")
        dump_file.write(",e_kdr,gbar_kdr")
        dump_file.write(",e_kdr_fs,gbar_kdr_fs")
        dump_file.write(",e_km,gbar_km")
        dump_file.write(",e_naf,gbar_naf")
        dump_file.write(",e_naf2,gbar_naf2")
        dump_file.write(",e_naf_tcr,gbar_naf_tcr")
        dump_file.write(",e_nap,gbar_nap")
        dump_file.write(",e_napf,gbar_napf")
        dump_file.write(",e_napf_spinstell,gbar_napf_spinstell")
        dump_file.write(",e_napf_tcr,gbar_napf_tcr")
        dump_file.write("\n")
        for comp in self.comp[1:]:
            dump_file.write('%s,%g,%g,%g,%g,%g,%g,%g,%g' % (comp.name, comp.length, comp.diameter, comp.sarea(), comp.xarea(), comp.Em, comp.Cm, comp.Rm, comp.Ra))
            for channel_name in config.channel_name_list:
                path = comp.path + '/' + channel_name
                if channel_name == 'CaPool':
                    if config.context.exists(path):
                        ca_pool = CaPool(path)
                        dump_file.write(',%g,%g' % (ca_pool.tau, ca_pool.B))
                    else:
                        dump_file.write(',0.0,0.0')
                    continue
                channel = None
                for channel in (ch for ch in comp.children() if ch.path().endswith(channel_name)):
                    break
                if channel is None:
                    dump_file.write(',0.0,0.0')
                else:
                    channel = moose.HHChannel(channel)
                    dump_file.write(',%g,%g' % (channel.Ek, channel.Gbar))
            dump_file.write('\n')



# 
# cell.py ends here
