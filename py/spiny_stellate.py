# spiny_stellate.py --- 
# 
# Filename: spiny_stellate.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Dec 19 14:26:52 2008 (+0530)
# Version: 
# Last-Updated: Sat Dec 20 16:38:55 2008 (+0530)
#           By: subhasis ray
#     Update #: 53
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

# Code:

from trb_tests import PyMooseTestContainer
from trb_channels import *

class SpinyStellate(PyMooseTestContainer):
    """Implementation of Spiny Stellate cell"""
    def __init__(self, *args):
        PyMooseTestContainer.__init__(self, *args)
        self.lib = moose.Neutral("/library")
        Globals.CONTEXT.setCwe(self.lib.id)
        
        # create the channel prototypes under /library
        self.naf_proto = NaF2("naf2")
        self.nap_proto = NaPSS("nap_ss")
        self.kdr_proto = KDRFS("kdr_fs")
        self.ka_proto = KA("ka")
        self.k2_proto = K2("k2")
        self.km_proto = KM("km")
        self.ar_proto = AR("ar")
        self.ca_conc_proto = CaConc("ca_conc")
        self.cal_proto = CaL("cal")
        self.cat_proto = CaT("cat")
        self.kc_proto = KCF("kc_fast")
        self.kahp_proto = KAHPSlow("kahp_slow")
        Globals.CONTEXT.setCwe("/")
        
        # Now read the cell
        Globals.CONTEXT.readCell("moose_spinystellate.p", self.model.path)
        
# 
# spiny_stellate.py ends here
