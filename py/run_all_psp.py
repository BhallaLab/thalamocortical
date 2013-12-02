#!/usr/bin/env python
"""Run all combinations of precell, postcell, synapse type for
plottting all the PSP"""

import os
import sys
from subprocess import call
from trbnetdata import TraubFullNetData
netdata = TraubFullNetData()
if __name__ == '__main__':
    # test_tcr_spinstell_ampa()
    # test_tcr_ss_spiking(0.0)
    for pretype in netdata.celltype:
        for posttype in netdata.celltype:
            for chantype in ['ampa', 'nmda', 'gaba']:
                call(['python', 'plot_psp.py', pretype, posttype, chantype])
