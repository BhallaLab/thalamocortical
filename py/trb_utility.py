# trb_utility.py --- 
# 
# Filename: trb_utility.py
# Description: utility functions
# Author: Subhasis Ray
# Maintainer: 
# Created: Sat Nov 29 04:41:16 2008 (+0530)
# Version: 
# Last-Updated: Fri Dec 19 11:14:51 2008 (+0530)
#           By: subhasis ray
#     Update #: 171
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

import logging

import moose

from trb_globals import *


def almost_equal(left, right, epsilon=0.001):
    """check if two floats are almost equal"""
    if left == right:
        return True
    if abs(left) > abs(right):
        return (1 - right / left) < epsilon
    else:
        return ( 1 - left / right) < epsilon
#!almost_equal


def list_subtree(root, node_list=[]):
    """Creates a list of all the objects in the subtree under root"""
    for child in root.children():
        dummy_obj = moose.Neutral(child)
        obj = eval("moose." + dummy_obj.className)(dummy_obj.id)
        node_list.append(obj)
        list_subtree(obj, node_list)
    return node_list
#! list_subtree


def compare_compartment(left, right):
    """Compare if two compartments have same field values"""
#     return almost_equal(left.Em, right.Em) and \
#         almost_equal(left.Rm, right.Rm) and \
#         almost_equal(left.Cm, right.Cm) and \
#         almost_equal(left.Ra, right.Ra) and \
#         almost_equal(left.initVm, right.initVm)
    result = almost_equal(left.Em, right.Em)
    if not result:
        logging.debug(left.path + ".Em = " + str(left.Em) + " <> " + right.path + ".Em = " + str(right.Em))
        return result
    result = almost_equal(left.Rm, right.Rm)
    if not result:
        logging.debug(left.path + ".Rm = " + str(left.Rm) + " <> " + right.path + ".Rm = " + str(right.Rm))
        return result

    result = almost_equal(left.Cm, right.Cm)
    if not result:
        logging.debug(left.path + ".Cm = " + str(left.Cm) + " <> " + right.path + ".Cm = " + str(right.Cm))
        return result
    result = almost_equal(left.Ra, right.Ra)
    if not result:
        logging.debug(left.path + ".Ra = " + str(left.Ra) + " <> " + right.path + ".Ra = " + str(right.Ra))
        return result

    result = almost_equal(left.initVm, right.initVm)
    if not result:
        logging.debug(left.path + ".initVm = " + str(left.initVm) + " <> " + right.path + ".initVm = " + str(right.initVm))
        return result

    return True
# ! compare_compartments
    

def compare_channel(left, right):
    """Compare two channels on same field values"""
#     result = almost_equal(left.Ek, right.Ek) and \
#         almost_equal(left.Gbar, right.Gbar) and \
#         almost_equal(left.Xpower, right.Xpower) and \
#         almost_equal(left.Ypower, right.Ypower) and \
#         almost_equal(left.Zpower, right.Zpower) and \
#         almost_equal(left.instant, right.instant) and \
#         almost_equal(left.Gk, right.Gk) and \
#         almost_equal(left.Ik, right.Ik)
    result = almost_equal(left.Ek, right.Ek)
    if not result:
        logging.debug(left.path + ".Ek = " + str(left.Ek) + " <> " + right.path + ".Ek = " + str(right.Ek))
        return result
    
    result = almost_equal(left.Gbar, right.Gbar)
    if not result:
        logging.debug(left.path + ".Gbar = " + str(left.Gbar) + " <> " + right.path + ".Gbar = " + str(right.Gbar))
        return result

    result = almost_equal(left.Xpower, right.Xpower)

    if not result:
        logging.debug(left.path + ".Xpower = " + str(left.Xpower) + " <> " + right.path + ".Xpower = " + str(right.Xpower))
        return result

    result = almost_equal(left.Ypower, right.Ypower)
    if not result:
        logging.debug(left.path + ".Ypower = " + str(left.Ypower) + " <> " + right.path + ".Ypower = " + str(right.Ypower))
        return result

    result = almost_equal(left.Zpower, right.Zpower)
    if not result:
        logging.debug(left.path + ".Zpower = " + str(left.Zpower) + " <> " + right.path + ".Zpower = " + str(right.Zpower))
        return result

    result = almost_equal(left.instant, right.instant)
    if not result:
        logging.debug(left.path + ".instant = " + str(left.instant) + " <> " + right.path + ".instant = " + str(right.instant))
        return result

    result = almost_equal(left.Gk, right.Gk)
    if not result:
        logging.debug(left.path + ".Gk = " + str(left.Gk) + " <> " + right.path + ".Gk = " + str(right.Gk))
        return result

    result = almost_equal(left.Ik, right.Ik)
    if not result:
        logging.debug(left.path + ".Ik = " + str(left.Ik) + " <> " + right.path + ".Ik = " + str(right.Ik))
        return result

    return True

# The following is unnecessary when we are comparing subtrees
#     if left.Xpower > 0:
#         result = result and \
#             compare_tables(moose.Table(left.path + "/xGate/A"), moose.Table(right.path + "/xGate/A")) and \
#             compare_tables(moose.Table(left.path + "/xGate/B"), moose.Table(right.path + "/xGate/B"))
#     if left.Ypower > 0:
#         result = result and \
#             compare_tables(moose.Table(left.path + "/yGate/A"), moose.Table(right.path + "/yGate/A")) and \
#             compare_tables(moose.Table(left.path + "/yGate/B"), moose.Table(right.path + "/yGate/B"))
#     return result
#! compare_channels

def compare_interpol(left, right):
    """Compare two interpol tables"""
    return left.mode == right.mode and \
        left.calcMode == right.calcMode and \
        left.xdivs == right.xdivs and \
        almost_equal(left.xmin, right.xmin) and \
        almost_equal(left.xmax, right.xmax) and \
        almost_equal(left.dx, right.dx) and \
        almost_equal(left.sy, right.sy) and \
        reduce( lambda x, y: x and y, map(almost_equal, left, right) )

def compare_table(left, right):
    """Compare two tables on same field values"""
    return left.stepMode == right.stepMode and \
        almost_equal(left.stepSize, right.stepSize) and \
        almost_equal(left.threshold, right.threshold) and \
        compare_interpol(left, right)

#! compare_tables


def compare_neutral(left, right):
    """Do nothing - no real parameters in neutral"""
    return True


def compare_element(left, right):
    """Compare element on the left to that on the right"""
    dummy_left = moose.Neutral(left.id)
    l_cls = dummy_left.className
    dummy_left = eval("moose." + l_cls)(left.id)

    dummy_right = moose.Neutral(right.id)
    r_cls = right.className
    dummy_right = eval("moose." + r_cls)(right.id) 
    
    if l_cls != r_cls:
        logging.debug("Objects of different class: " + left.className + ":" + left.path + " <> " + right.className + ":" + right.path)
        return False
    method_ = compare_neutral
    if l_cls == "Compartment":
        method_ = compare_compartment
    elif l_cls == "HHChannel":
        method_ = compare_channel
    elif l_cls == "Interpol":
        method_ = compare_interpol
    elif l_cls == "Table":
        method_ = compare_table

    return method_(dummy_left, dummy_right)
# ! compare_element


def compare_subtree(left, right):
    """Compare two model subtrees rooted at left and right.

    It returns a 3-tuple containing the result, the element on the
    left and that on the right which violate equality.

    The limitation is that the names must match of every element in
    left with every element in right subtree."""
    result = True
    result = result and compare_element(left, right)
    if result is False:
        return (result, left, right)
    left_children = []
    right_children = []
    for left_child in left.children():
        left_children.append(left_child)
    for right_child in right.children():
        right_children.append(right_child)
    if len(left_children) != len(right_children):
        logging.debug("Number of children not matching: " + left.path + ": " + str(len(left_children)) + " vs " + right.path + ": " + str(len(right_children)))
        return (False, left, right)
    # I am not sure if the children() method will return a
    # deterministic ordering of children independent of their creation
    # time (Id). Hence better be sure by sorting on the path
    left_children.sort(key=lambda x: x.path())
    right_children.sort(key=lambda x: x.path())
    for child in left_children:
        print child.path(),
    print ""
    for child in right_children:
        print child.path(),
    print ""
    for left_child, right_child in zip(left_children, right_children):
        (new_result, lft, rht) = compare_subtree(moose.Neutral(left_child), moose.Neutral(right_child))
#         print lft.path, rht.path, result
        result = result and new_result
        if result is False:
            return (result, lft, rht)
    return (result, None, None)
#! compare_subtree

def insert_injection(comp_path, delay, width, level):
    """Insert a pulsegen into a compartment as a source of current
    injection."""
    pulsegen = moose.PulseGen(comp_path + "/injection")
    pulsegen.firstDelay = delay
    pulsegen.firstWidth = width
    pulsegen.firstLevel = level
    pulsegen.connect("outputSrc", comp, "injectMsg")
    return pulsegen


def init_changates(channel, xPower, yPower = 0):
    """Initialize xGate and yGate and the tables inside them as required."""
    gates = []
    if xPower != 0:
        channel.Xpower = xPower
        xGate = moose.HHGate(channel.path+"/xGate")
        xGate.A.xmin = Globals.VMIN
        xGate.A.xmax = Globals.VMAX
        xGate.A.xdivs = Globals.NDIVS
        xGate.B.xmin = Globals.VMIN
        xGate.B.xmax = Globals.VMAX
        xGate.B.xdivs = Globals.NDIVS
        gates.append(xGate)

    if yPower != 0:
        channel.Ypower = yPower
        yGate = moose.HHGate(channel.path+"/yGate")
        yGate.A.xmin = Globals.VMIN
        yGate.A.xmax = Globals.VMAX
        yGate.A.xdivs = Globals.NDIVS
        yGate.B.xmin = Globals.VMIN
        yGate.B.xmax = Globals.VMAX
        yGate.B.xdivs = Globals.NDIVS
        gates.append(yGate)
    return gates

def create_channel_output(channel, data_container):
    """Create recording tables for the specified channel"""
    parent = moose.Neutral(channel.parent)
    ik_table = moose.Table(parent.name + "_" + channel.name + "_Ik", data_container)
    ik_table.stepMode = 3
    ik_table.xmin = 0.0
    ik_table.xmax = 1.0
    ik_table.xdivs = int(Globals.plotsteps)
    ik_table.connect("inputRequest", channel, "Ik")
    gk_table = moose.Table(parent.name + "_" + channel.name + "_Gk", data_container)
    gk_table.stepMode = 3
    gk_table.xmin = 0.0
    gk_table.xmax = 1.0
    gk_table.xdivs = int(Globals.plotsteps)
    gk_table.connect("inputRequest", channel, "Gk")
#!create_channel_output    
    
                           
# 
# trb_utility.py ends here
