# spinystellate.py --- 
# 
# Filename: spinystellate.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Tue Sep 29 11:43:22 2009 (+0530)
# Version: 
# Last-Updated: Tue Sep 29 11:59:21 2009 (+0530)
#           By: subhasis ray
#     Update #: 8
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

from datetime import datetime
import moose
import config
from cell import *
from capool import CaPool

class SpinyStellate(TraubCell):
    prototype = TraubCell.read_proto("SpinyStellate.p", "SpinyStellate")
    def __init__(self, *args):
	TraubCell.__init__(self, *args)

    def _topology(self):
	self.presyn = 57

    def _setup_passive(self):
	for comp in self.comp[1:]:
	    comp.initVm = -65e-3

    def _setup_channels(self):


# 
# spinystellate.py ends here
