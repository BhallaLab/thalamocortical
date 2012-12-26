# trbutil.py --- 
# 
# Filename: trbutil.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Fri Jun  5 13:59:40 2009 (+0530)
# Version: 
# Last-Updated: Wed Dec 26 09:41:55 2012 (+0530)
#           By: subha
#     Update #: 17
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

import smtplib
from subprocess import call
import numpy
import config
import gzip
import numpy
from scipy import signal

def almost_equal(left, right, epsilon=1e-6):
    """check if two floats are almost equal"""
    if left == right:
        return True
    if abs(left) > abs(right):
        return (1 - right / left) < epsilon
    else:
        return ( 1 - left / right) < epsilon
#!almost_equal

def read_nrn_data(filename, hoc_script=None):
    data = None
    filepath = '../nrn/data/' + filename
    try:
        data = numpy.loadtxt(filepath)
    except IOError:
        try:
            data = numpy.loadtxt(filepath + '.gz')
        except IOError:
            call([config.neuron_bin, hoc_script], cwd='../nrn')
            data = numpy.loadtxt(filepath)
    return data

def do_plot(class_name, mus_t, mus_ca, mus_vm, nrn_t=None, nrn_ca=None, nrn_vm=None):
    '''Plot the membrane potential over time in both moose and neuron.'''
    if has_pylab:
        if nrn_vm is None or len(nrn_vm) is 0:
            nrn_t = pylab.array()
            nrn_vm = pylab.array() 
            nrn_ca  = pylab.array()
        config.pylab.subplot(211)
        config.pylab.plot(nrn_t, nrn_vm, 'y-', label='NEURON')
        config.pylab.plot(mus_t, mus_vm, 'g-.', label='MOOSE')
        config.pylab.title('Vm in presynaptic compartment of %s' % class_name)
        config.pylab.legend()
        config.pylab.subplot(212)
        config.pylab.plot(nrn_t, nrn_ca, 'r-', label='NEURON')
        config.pylab.plot(mus_t, mus_ca, 'b-.', label='MOOSE')
        config.pylab.title('[Ca2+] in soma of %s' % class_name)
        config.pylab.legend()
        config.pylab.show()

def send_email(recipient, sender, password,
               server='smtp.gmail.com', 
               port=587, subject='no subject', body='empty message'):
    """Send an email to specified `recipient` using `sender`
    account on `server`"""
    smtp_obj = smtplib.SMTP(server, port)
    smtp_obj.starttls()
    smtp_obj.login(sender, password)
    headers = ["from: " + sender,
           "subject: " + subject,
               "to: " + recipient,
               "mime-version: 1.0",
               "content-type: text/html"]
    headers = "\r\n".join(headers)
    smtp_obj.sendmail(sender, recipient, headers + "\r\n\r\n" + body)
    smtp_obj.quit()

if __name__ == '__main__':
    read_nrn_data('asa')

# 
# trbutil.py ends here
