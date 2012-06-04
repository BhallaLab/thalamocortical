# test_naf.py --- 
# 
# Filename: test_naf.py
# Description: 
# Author: 
# Maintainer: 
# Created: Sat May 26 12:02:48 2012 (+0530)
# Version: 
# Last-Updated: Thu May 31 15:41:51 2012 (+0530)
#           By: subha
#     Update #: 78
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

import uuid
import unittest
import numpy as np
import moose
import nachans
import kchans

from testutils import *
library = moose.Neutral('/library')

class ChannelTestBase(unittest.TestCase):    
    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        self.uuid = uuid.uuid4().int
        self.container = moose.Neutral('test%d' % (self.uuid))
        self.params = setup_single_compartment(self.container.path, self.get_protochan())
        moose.context.reset()
        moose.context.step(0.35)
        vm_file = '%s_Vm.dat' % (self.get_protochan().name)
        gk_file = '%s_Gk.dat' % (self.get_protochan().name)
        self.params['Vm'].dumpFile(vm_file)
        print 'Saved Vm in', vm_file
        self.params['Gk'].dumpFile(gk_file)
        print 'Saved Gk in', gk_file
        
    def get_protochan(self):
        raise NotImplementedError('Subclasses must implement this method')
    
        
class TestNaF(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaF('NaF', library)
        
    def testNaF(self):
        pass


class TestNaF2(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaF2('NaF2', library)
    
    def testNaF2(self):
        pass

    
class TestNaF2_nRT(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaF2_nRT('NaF2_nRT', library)
    
    def testNaF2_nRT(self):
        pass

    
class TestNaP(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaP('NaP', library)

    def testNaP(self):
        pass


class TestNaPF(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaPF('NaPF', library)

    def testNaPF(self):
        pass


class TestNaPF_SS(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaPF_SS('NaPF_SS', library)

    def testNaPF_SS(self):
        pass

    
class TestNaPF_TCR(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaPF_TCR('NaPF_TCR', library)

    def testNaPF_TCR(self):
        pass


class TestNaF_TCR(ChannelTestBase):
    def get_protochan(self):
        return nachans.NaF_TCR('NaF_TCR', library)

    def testNaF_TCR(self):
        pass

class TestK2(ChannelTestBase):
    def get_protochan(self):
        return kchans.K2('K2', library)

    def testK2(self):
        pass
    
    
if __name__ == '__main__':
    unittest.main()
    
# 
# test_naf.py ends here
