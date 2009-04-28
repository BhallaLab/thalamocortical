# compartment.py --- 
# 
# Filename: compartment.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Apr 24 10:01:45 2009 (+0530)
# Version: 
# Last-Updated: Tue Apr 28 20:58:35 2009 (+0530)
#           By: subhasis ray
#     Update #: 41
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
import moose
import config
from math import pi
from cachans import *
from kchans import *
from nachans import *
from capool import *
from archan import AR

class MyCompartment(moose.Compartment):
    def __init__(self, *args):
        moose.Compartment.__init__(self, *args)
        self.channels = []
        self._xarea = None
        self._sarea = None

    def setSpecificRm(self, RM):
        self.Rm = RM / self.sarea()

    def setSpecificRa(self, RA):
        self.Ra = RA * self.length / self.xarea()

    def setSpecificCm(self, CM):
        self.Cm = CM * self.sarea()

    def xarea(self):
        if self._xarea is None:
            self._xarea = pi * self.diameter * self.diameter
        return self._xarea

    def sarea(self):
        if self._sarea is None:
            self._sarea = pi * self.length * self.diameter
        return self._sarea

    def insertChannel(self, channel, specificGbar=None, Ek=None):
        """Insert a channel setting its gbar as membrane_area *
        specificGbar and reversal potential to Ek.
        
        This method expects either a valid channel class name or an
        existing channel object. If specificGbar is given, the Gbar is
        set to specificGbar * surface-area of the compartment. If Ek
        is given, the channel's Ek is set to this value.
        """
        if type(channel) is type(''): # if it is a class name, create the channel as a child with the same name as the class name
            
            chan_class = eval(channel)
            chan = chan_class(channel, self)
        elif type(channel) is moose.HHChannel:
            chan = channel
        else:
            print "ERROR: unknown object passed as channel: ", channel
        if specificGbar is not None:
            chan.Gbar = specificGbar * self.sarea()
        if Ek is not None:
            chan.Ek = Ek
        self.channels.append(chan)
        self.connect("channel", chan, "channel")
        return chan

    def insertCaPool(self, phi, tau):
        """Insert a Ca+2 pool and connect it to the relevant channels.

        phi is the amount of Ca2+ in unit area. 

        NOTE that this function should be called only after all
        channels (Ca and Ca dependent K channels) have been
        initialized. You can call this function multiple times without
        harm, there is safeguard against multiple connections in
        CaPool class."""
        self.ca_pool = CaPool('CaPool', self)
        self.ca_pool.B = phi / self.sarea()
        self.ca_pool.tau = tau
        ca_channels = [ channel for channel in self.channels \
                            if isinstance(channel, CaL) ]
        self.ca_pool.connectCaChannels(ca_channels)
        kca_channels = [ channel for channel in self.channels \
                            if isinstance(channel, KCaChannel) ]
        self.ca_pool.connectDepChannels(kca_channels)

    def insertRecorder(self, field_name, data_container):
        """Creates a table for recording a field under data_container"""
        table = moose.Table(field_name, data_container)
        table.stepMode = 3
        self.connect(field_name, table, "inputRequest")
        return table

    def insertPulseGen(self, name, parent, baseLevel=0.0, firstLevel=1e-10, firstDelay=20e-3, firstWidth=20e-3, secondLevel=0.0, secondDelay=1e10, secondWidth=0.0):
        self.pulsegen = moose.PulseGen(name, parent)
        self.pulsegen.baseLevel = baseLevel
        self.pulsegen.firstLevel = firstLevel
        self.pulsegen.firstDelay = firstDelay
        self.pulsegen.firstWidth = firstWidth
        self.pulsegen.secondLevel = secondLevel
        self.pulsegen.secondDelay = secondDelay
        self.pulsegen.secondWidth = secondWidth
        self.pulsegen.connect('outputSrc', self, 'injectMsg')
        return self.pulsegen





        
# 
# compartment.py ends here
