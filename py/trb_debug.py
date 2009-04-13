# trb_debug.py --- 
# 
# Filename: trb_debug.py
# Description: Equivaent to debug.g. 
#              Contains functions for logging fields of an objects.
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Nov 29 03:43:29 2008 (+0530)
# Version: 
# Last-Updated: Mon Dec  8 17:50:56 2008 (+0530)
#           By: subhasis ray
#     Update #: 58
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

import logging

import moose
from trb_globals import *

def log_fields(obj, field_list):
    """Log the fields to the log file"""
    for field in field_list:
        logging.info(obj.path + "." + field + " = " + str(getattr(obj, field)))
# !log_fields


def dump_tables(container):
    """Save contents of each table in a separate file.

    Container is the object within which the table objects are located.
    """
    for child in container.children():
        table = moose.Table(child)
        file_name = Globals.SIMULATOR + "_" + table.name + ".plot"
        table.dumpFile(file_name)

# A map of channel gate table paths and corresponding file name specializer
gate_table_map = {"xGate/A": "xa",
                  "xGate/B": "xb",
                  "yGate/A": "ya",
                  "yGate/B": "yb"}


def dump_channel_tables(channel):
    """Dump the channel-gate tables in files."""
    for key in gate_table_map.keys():
        path = channel.path + "/" + key
        file_path = Globals.SIMULATOR + "_" + channel.name + "_" + gate_table_map[key]+ ".plot"
        if channel.getContext().exists(path):
            print "dumped: ", path, "in", file_path
            moose.Table(path).dumpFile(file_path)
        else:
            print path, "does not exist"
#!save_tables

# 
# trb_debug.py ends here
