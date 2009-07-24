# cell.py --- 
# 
# Filename: cell.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Jul 24 10:04:47 2009 (+0530)
# Version: 
# Last-Updated: Fri Jul 24 14:07:02 2009 (+0530)
#           By: subhasis ray
#     Update #: 16
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
import config
import moose

class TraubCell(moose.Cell):
    def __init__(self, *args):
	moose.Cell.__init__(self, *args)
	self.comp = [None]
	self.num_comp = 0
        self.level = defaultdict(set)
	self.dendrites = set()
        self.presyn = 0

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
                if channel_name == 'CaPool':
                    if hasattr(comp, 'ca_pool'):
                        dump_file.write(',%g,%g' % (comp.ca_pool.tau, comp.ca_pool.B * comp.sarea()))
                    else:
                        dump_file.write(',0.0,0.0')
                    continue
                channel = None
                for channel in (ch for ch in comp.channels if ch.name == channel_name):
                    break
                if channel is None:
                    dump_file.write(',0.0,0.0')
                else:
                    dump_file.write(',%g,%g' % (channel.Ek, channel.Gbar))
            dump_file.write('\n')



# 
# cell.py ends here
