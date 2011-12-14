# single_cell_regression.py --- 
# 
# Filename: single_cell_regression.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Mon Dec 12 11:52:46 2011 (+0530)
# Version: 
# Last-Updated: Wed Dec 14 10:58:20 2011 (+0530)
#           By: Subhasis Ray
#     Update #: 14
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# Run a regression test on all single cell models.
# 
# 

# Change log:
# 
# 
# 

# Code:

import subprocess

singlecell_files = [
    'deepaxoaxonic.py',
    'deepbasket.py',
    'deepLTS.py',
    'nontuftedRS.py',
    'nRT.py',
    'spinystellate.py',
    'supaxoaxonic.py',
    'supbasket.py',
    'supLTS.py',
    'suppyrFRB.py',
    'suppyrRS.py',
    'tcr.py',
    'tuftedIB.py',
    'tuftedRS.py']
if __name__ == '__main__':
    for filename in singlecell_files:
        subprocess.call(['python', filename])


# 
# single_cell_regression.py ends here
