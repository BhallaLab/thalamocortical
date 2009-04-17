# trb_tests.py --- 
# 
# Filename: trb_tests.py
# Description: Contains unit test code for Traub 2005
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Dec  2 23:03:52 2008 (+0530)
# Version: 
# Last-Updated: Mon Apr 13 16:25:51 2009 (+0530)
#           By: subhasis ray
#     Update #: 402
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

from math import *
import unittest
import logging
import pylab

import sys

sys.path.append("/home/subha/src/moose/pymoose")

import moose

from trb_globals import Globals
from trb_debug import *
from trb_utility import *
from trb_channels import *

TEST_RA = 2.5
TEST_RM = 5.0
TEST_CM = 0.009
TEST_EM = -68e-3
TEST_EK = -100e-3
TEST_DIA = 15e-6
TEST_LEN = 20e-6

TEST_GNAF_DENS = 1500.0
TEST_GKDR_DENS = 1000.0
TEST_GNAP_DENS = 0.001 * TEST_GNAF_DENS
TEST_GCAT_DENS = 1.0
TEST_GCAL_DENS = 5.0
TEST_GKA_DENS = 300.0
TEST_GKC_DENS = 100.0
TEST_GKM_DENS = 37.5
TEST_GK2_DENS = 1.0
TEST_GKAHP_DENS = 1.0
TEST_GAR_DENS = 2.5

PI = 3.141592

TEST_SURFACE = PI * TEST_DIA * TEST_LEN
TEST_XAREA = PI * TEST_DIA * TEST_DIA / 4.0

PRE_TIME = 0.05
PULSE_WIDTH = 0.05
POST_TIME = 0.05

class TestCompartment(moose.Compartment):
    """A compartment with default field values."""
    def __init__(self, *args):
        moose.Compartment.__init__(self, *args)
        self.diameter = TEST_DIA
        self.length = TEST_LEN
        self.Em = TEST_EM
        self.initVm = TEST_EM
        self.Ra = TEST_RA * TEST_LEN / TEST_XAREA
        self.Rm = TEST_RM / TEST_SURFACE
        self.Cm = TEST_CM * TEST_SURFACE
    #! __init__
#! class TestCompartment


class PyMooseTestContainer(moose.Neutral):
    """A container for tests. 

    It has two children: model and data.  model contains the model
    tree. data contains the tables for recording time series data."""
    def __init__(self, *args):
        moose.Neutral.__init__(self, *args)
        self.model = moose.Neutral("model", self)
        self.data = moose.Neutral("data", self)

    def schedule(self):
        """Set the clocks and schedule objects to clocks"""
        Globals.CONTEXT.setClock(0, Globals.SIMDT, 0)
        Globals.CONTEXT.setClock(1, Globals.SIMDT, 1)        
        Globals.CONTEXT.setClock(2, Globals.PLOTDT, 0)
        Globals.CONTEXT.useClock(0, self.model.path + "/##")
        Globals.CONTEXT.useClock(2, self.data.path + "/#")
        
    def run_simulation(self):
        """Reset and step. Finally dump the data"""
        Globals.CONTEXT.reset()
        Globals.CONTEXT.step(Globals.simtime)

    def dump_data(self):
        """Save all the tables in ./data to files"""
        dump_tables(self.data)

# ! PyMooseTestContainer


class SingleCompPassive(PyMooseTestContainer):
    """Creates full model with a passive single compartment and the recording tables"""
    def __init__(self, *args):
        PyMooseTestContainer.__init__(self, *args)
        self.comp = TestCompartment("comp", self.model)
        self.vm_table = moose.Table("Vm", self.data)
        self.vm_table.stepMode = 3
        self.vm_table.connect("inputRequest", self.comp, "Vm")

        self.pulsegen = moose.PulseGen("inject", self.model)
        self.pulsegen.firstDelay = PRE_TIME
        self.pulsegen.firstWidth = PULSE_WIDTH
        self.pulsegen.firstLevel = Globals.INJECTION
        self.pulsegen.connect("outputSrc", self.comp, "injectMsg")

        self.inject_table = moose.Table("Inject", self.data)
        self.inject_table.stepMode = 3
        self.inject_table.connect("inputRequest", self.pulsegen, "output")
    #! __init__
# !SingleCompPassive


class SingleCompInjectTestCase(unittest.TestCase):
    """Base class for test cases with a single compartment with current injection"""
    def __init__(self, *args):
        unittest.TestCase.__init__(self, *args)
        Globals.simtime = PRE_TIME + PULSE_WIDTH + POST_TIME
        Globals.plotsteps = Globals.simtime / Globals.PLOTDT
        Globals.simsteps = Globals.simsteps / Globals.SIMDT
        self._test_obj = None

    def get_test_object(self, root="pymoose_single_comp_passive"):
        if self._test_obj is None:
            self._test_obj = SingleCompPassive(root)
        return self._test_obj
# ! SingleCompInjectTestCase         


class SingleCompPassiveTestCase(SingleCompInjectTestCase):
    """Test a single passive compartment"""
    def __init__(self, *args):
        SingleCompInjectTestCase.__init__(self, *args)

    #! __init__
    
    def testSingleCompPassive(self):
        """Run the simulation and dump the data"""
        self.get_test_object("pymoose_single_comp_passive").schedule()
        self.get_test_object().run_simulation()
        self.get_test_object().dump_data()
    # !testSingleCompPassive
# ! SingleCompPassiveTestCase


class MOOSEPySingleCompPassiveTestCase(unittest.TestCase):
    """Create the MOOSE and PyMOOSE models and compare them"""
    def __init__(self, *args):
        unittest.TestCase.__init__(self, *args)
        self.pymoose_root = "pymoose_single_comp_passive"
        self.moose_root = "moose_single_comp_passive"
        self.pymoose_model = SingleCompPassive(self.pymoose_root).model
        Globals.CONTEXT.loadG("trb_tests.g")
        Globals.CONTEXT.runG("create_test_container " + self.moose_root)
        Globals.CONTEXT.runG("setup_single_compartment_passive " + self.moose_root)
        self.moose_model = moose.Neutral(self.moose_root + "/model")
    # ! __init__

    def testSubtreeEquality(self):
        """Compare objects field by field"""
        self.assert_(compare_subtree(self.pymoose_model, self.moose_model))
#! class MOOSEPySingleCompPassiveTestCase


class SingleCompMultiChannel(SingleCompPassive):
    """A single compartment with multiple channels"""
    def __init__(self, *args):
        SingleCompPassive.__init__(self, *args)

        self.naf = NaF2("naf", self.comp)
        self.naf.Ek = Globals.E_NA
        self.naf.Gbar = TEST_GNAF_DENS * TEST_SURFACE
        self.naf.connect("channel", self.comp, "channel")
        dump_channel_tables(self.naf)
        create_channel_output(self.naf, self.data)

        self.nap = NaPSS("nap", self.comp)
        self.nap.Ek = Globals.E_NA
        self.nap.Gbar = TEST_GNAP_DENS * TEST_SURFACE
        self.nap.connect("channel", self.comp, "channel")
        dump_channel_tables(self.nap)
        create_channel_output(self.nap, self.data)

        self.kdr = KDRFS("kdr", self.comp)
        self.kdr.Ek = Globals.E_K_FS
        self.kdr.Gbar = TEST_GKDR_DENS * TEST_SURFACE
        self.kdr.connect("channel", self.comp, "channel")
        dump_channel_tables(self.kdr)
        log_fields(self.kdr, Globals.HHCHANNEL_FIELDS)
        create_channel_output(self.kdr, self.data)

        self.ka = KA("ka", self.comp)
        self.ka.Ek = Globals.E_K_FS
        self.ka.Gbar = TEST_GKA_DENS * TEST_SURFACE
        self.ka.connect("channel", self.comp, "channel")
        dump_channel_tables(self.ka)
        create_channel_output(self.ka, self.data)

        self.k2 = K2("k2", self.comp)
        self.k2.Ek = Globals.E_K_FS
        self.k2.Gbar = TEST_GK2_DENS * TEST_SURFACE
        self.k2.connect("channel", self.comp, "channel")
        dump_channel_tables(self.k2)
        create_channel_output(self.k2, self.data)

        self.km = KM("km", self.comp)
        self.km.Ek = Globals.E_K_FS
        self.km.Gbar = TEST_GKM_DENS * TEST_SURFACE
        self.km.connect("channel", self.comp, "channel")
        dump_channel_tables(self.km)
        create_channel_output(self.km, self.data)

        self.ar = AR("ar", self.comp)
        self.ar.Ek = Globals.E_AR
        self.ar.Gbar = TEST_GAR_DENS * TEST_SURFACE
        self.ar.connect("channel", self.comp, "channel")
        dump_channel_tables(self.ar)
        create_channel_output(self.ar, self.data)

        self.ca_conc = moose.CaConc("ca_conc", self.comp)
        self.ca_conc.B = 1.0 / (2 * 96485 * PI * 0.2e-6 * self.comp.diameter * self.comp.length)
        self.ca_conc.tau = 50e-3

        self.cal = CaL("cal", self.comp)
        self.cal.Ek = Globals.E_CA
        self.cal.Gbar = TEST_GCAL_DENS * TEST_SURFACE
        self.cal.connect("channel", self.comp, "channel")
        self.cal.connect("IkSrc", self.ca_conc, "current")
        dump_channel_tables(self.cal)
        create_channel_output(self.cal, self.data)

        self.cat = CaTA("cat", self.comp)
        self.cat.Ek = Globals.E_CA
        self.cat.Gbar = TEST_GCAT_DENS * TEST_SURFACE
        self.cat.connect("channel", self.comp, "channel")
        self.cat.connect("IkSrc", self.ca_conc, "current")
        dump_channel_tables(self.cat)
        create_channel_output(self.cat, self.data)

        self.kc = KC("kc", self.comp)
        self.kc.Ek = Globals.E_K_FS
        self.kc.Gbar = TEST_GKC_DENS * TEST_SURFACE
        self.kc.connect("channel", self.comp, "channel")
        dump_channel_tables(self.kc)
        create_channel_output(self.kc, self.data)
        self.ca_conc.connect("concSrc", self.kc, "concen")

        self.kahp = KAHPSlow("kahp", self.comp)
        self.kahp.Ek = Globals.E_K_FS
        self.kahp.Gbar = TEST_GKAHP_DENS * TEST_SURFACE
        self.kahp.connect("channel", self.comp, "channel")
        dump_channel_tables(self.kahp)
        create_channel_output(self.kahp, self.data)
        self.ca_conc.connect("concSrc", self.kahp, "concen")
    # !__init__
#! class SingleCompMultiChannel


class SingleCompMultiChannelTestCase(SingleCompInjectTestCase):
    """Test a single compartment with multiple channels"""
    def __init__(self, *args):
        SingleCompInjectTestCase.__init__(self, *args)
        
    def get_test_object(self, root="pymoose_single_comp_multichan"):
        """Overrides the method in SingleCompInjectTestCase.
        
        creates the _test_obj member variable (if not present) as an instance of SingleCompMultiChannel.
        returns the _test_obj member variable.
        """
        if self._test_obj is None:
            self._test_obj = SingleCompMultiChannel(root)
        return self._test_obj

    def testSingleCompMultiChannel(self):
        """Run the simulation and dump the recorded time series in files."""
        self.get_test_object().schedule()
        self.get_test_object().run_simulation()
        self.get_test_object().dump_data()
    # ! testSingleCompMultiChannel
# ! SingleCompMultiChannelTestCase


class MOOSEPySingleCompMultiChannelTestCase(unittest.TestCase):
    """Create the MOOSE and PyMOOSE models and compare them"""
    def __init__(self, *args):
        unittest.TestCase.__init__(self, *args)
        self.moose_root = "moose_single_comp_multichan"
        self.pymoose_root = "pymoose_single_comp_multichan"
        self.pymoose_model = SingleCompMultiChannel(self.pymoose_root).model
        log_fields(moose.HHChannel(self.pymoose_model.path + "/comp/kdr"), Globals.HHCHANNEL_FIELDS)
        Globals.CONTEXT.loadG("trb_tests.g")
        Globals.CONTEXT.runG("create_test_container " + self.moose_root)
        Globals.CONTEXT.runG("setup_single_compartment_multichannel " + self.moose_root)
        self.moose_model = moose.Neutral(self.moose_root + "/model")
        log_fields(moose.HHChannel(self.pymoose_model.path + "/comp/kdr"), Globals.HHCHANNEL_FIELDS)
#         l1 = pylab.plot(moose.Interpol(self.pymoose_model.path + "/comp/naf/yGate/A"))
#         l2 = pylab.plot(moose.Interpol(self.moose_model.path + "/comp/naf/yGate/A"))
#         pylab.figlegend((l1, l2), ("pymoose: NaF Y_A", "moose: NaF Y_A"))
#         pylab.show()
#         l1 = pylab.plot(moose.Interpol(self.pymoose_model.path + "/comp/naf/yGate/B"))
#         l2 = pylab.plot(moose.Interpol(self.moose_model.path + "/comp/naf/yGate/B"))
#         pylab.figlegend((l1, l2), ("pymoose: NaF Y_B", "moose: NaF Y_B"))
#         pylab.show()
    # ! __init__

    def testSubtreeEquality(self):
        """Compare objects field by field"""
        (result, left, right) = compare_subtree(self.pymoose_model, self.moose_model)
        l_path = "None"
        r_path = "None"
        if left is not None:
            l_path = left.path
        if right is not None:
            r_path = right.path
        self.assert_(result, "Not matching: " + l_path + " and " + r_path)

#! class MOOSEPySingleCompMultiChannelTestCase



if __name__ == "__main__":
# Section 1: run the multichannel single compartment model
#     suite = unittest.TestLoader().loadTestsFromTestCase(SingleCompMultiChannelTestCase)
#     unittest.TextTestRunner(verbosity=2).run(suite)

# Section  2: test equivalence of the MOOSE model and the PyMOOSE model
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(filename)-10s: %(lineno)-6d: %(levelname)-8s %(message)s",filename="pymoose.log",filemode="w")
    suite = unittest.TestLoader().loadTestsFromTestCase(MOOSEPySingleCompMultiChannelTestCase)
    unittest.TextTestRunner(verbosity=2).run(suite)

# Run the single passive compartment model
#     testCase = MOOSEPySingleCompPassiveTestCase()
#     testCase.testParameters()

# 
# trb_tests.py ends here
