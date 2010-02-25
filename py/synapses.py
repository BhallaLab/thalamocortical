# synapses.py --- 
# 
# Filename: synapses.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Feb 25 15:22:11 2010 (+0530)
# Version: 
# Last-Updated: Fri Feb 26 01:33:54 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 39
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This file is for data on synapses of various
# kinds. Ideally it should be replaced by something in netCDF/HDF5 or
# some other data format. But given the small size of this data, I
# don't see much savings in that.
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

from collections import defaultdict

TAU_AMPA = defaultdict(dict)
TAU_NMDA = defaultdict(dict)
TAU_GABA = defaultdict(dict)
TAU_GABA_FAST = defaultdict(dict)
TAU_GABA_SLOW = defaultdict(dict)
G_AMPA = defaultdict(dict)
G_NMDA = defaultdict(dict)
G_GABA = defaultdict(dict)

# Synaptic conductance time constants. 
          tauAMPA_suppyrRS_to_suppyrRS=2.e0 
          tauNMDA['SupPyrRS']to_suppyrRS=130.5e0 
          tauAMPA['SupPyrRS']to_suppyrFRB=2.e0 
          tauNMDA['SupPyrRS']to_suppyrFRB=130.e0  
          tauAMPA['SupPyrRS']to_supbask  =.8e0   
          tauNMDA['SupPyrRS']to_supbask  =100.e0 
          tauAMPA['SupPyrRS']to_supaxax  =.8e0  
          tauNMDA['SupPyrRS']to_supaxax  =100.e0 
          tauAMPA['SupPyrRS']to_supLTS   =1.e0  
          tauNMDA['SupPyrRS']to_supLTS   =100.e0 
          tauAMPA['SupPyrRS']to_spinstell=2.e0   
          tauNMDA['SupPyrRS']to_spinstell=130.e0 
          tauAMPA['SupPyrRS']to_tuftIB   =2.e0   
          tauNMDA['SupPyrRS']to_tuftIB   =130.e0 
          tauAMPA['SupPyrRS']to_tuftRS   =2.e0   
          tauNMDA['SupPyrRS']to_tuftRS   =130.e0 
          tauAMPA['SupPyrRS']to_deepbask =.8e0   
          tauNMDA['SupPyrRS']to_deepbask =100.e0 
          tauAMPA['SupPyrRS']to_deepaxax =.8e0   
          tauNMDA['SupPyrRS']to_deepaxax =100.e0 
          tauAMPA['SupPyrRS']to_deepLTS  =1.e0   
          tauNMDA['SupPyrRS']to_deepLTS  =100.e0 
          tauAMPA['SupPyrRS']to_nontuftRS=2.e0   
          tauNMDA['SupPyrRS']to_nontuftRS=130.e0 

          tauAMPA['SupPyrFRB']to_suppyrRS=2.e0   
          tauNMDA['SupPyrFRB']to_suppyrRS=130.e0  
          tauAMPA['SupPyrFRB']to_suppyrFRB=2.e0   
          tauNMDA['SupPyrFRB']to_suppyrFRB=130.e0 
          tauAMPA['SupPyrFRB']to_supbask  =.8e0   
          tauNMDA['SupPyrFRB']to_supbask  =100.e0 
          tauAMPA['SupPyrFRB']to_supaxax  =.8e0  
          tauNMDA['SupPyrFRB']to_supaxax  =100.e0 
          tauAMPA['SupPyrFRB']to_supLTS   =1.e0  
          tauNMDA['SupPyrFRB']to_supLTS   =100.e0 
          tauAMPA['SupPyrFRB']to_spinstell=2.e0   
          tauNMDA['SupPyrFRB']to_spinstell=130.e0 
          tauAMPA['SupPyrFRB']to_tuftIB   =2.e0   
          tauNMDA['SupPyrFRB']to_tuftIB   =130.e0 
          tauAMPA['SupPyrFRB']to_tuftRS   =2.e0   
          tauNMDA['SupPyrFRB']to_tuftRS   =130.e0
          tauAMPA['SupPyrFRB']to_deepbask =.8e0   
          tauNMDA['SupPyrFRB']to_deepbask =100.e0 
          tauAMPA['SupPyrFRB']to_deepaxax =.8e0   
          tauNMDA['SupPyrFRB']to_deepaxax =100.e0 
          tauAMPA['SupPyrFRB']to_deepLTS  =1.e0   
          tauNMDA['SupPyrFRB']to_deepLTS  =100.e0 
          tauAMPA['SupPyrFRB']to_nontuftRS=2.e0   
          tauNMDA['SupPyrFRB']to_nontuftRS=130.e0 

          tauGABA['SupBasket']to_suppyrRS   =6.e0  
          tauGABA['SupBasket']to_suppyrFRB  =6.e0  
          tauGABA['SupBasket']to_supbask    =3.e0  
          tauGABA['SupBasket']to_supaxax    =3.e0  
          tauGABA['SupBasket']to_supLTS     =3.e0  
          tauGABA['SupBasket']to_spinstell  =6.e0  

          tauGABA['SupAxoaxonic']to_suppyrRS   =6.e0  
          tauGABA['SupAxoaxonic']to_suppyrFRB  =6.e0  
          tauGABA['SupAxoaxonic']to_spinstell  =6.e0  
          tauGABA['SupAxoaxonic']to_tuftIB     =6.e0  
          tauGABA['SupAxoaxonic']to_tuftRS     =6.e0  
          tauGABA['SupAxoaxonic']to_nontuftRS  =6.e0  

          tauGABA['SupLTS']to_suppyrRS    =20.e0 
          tauGABA['SupLTS']to_suppyrFRB   =20.e0 
          tauGABA['SupLTS']to_supbask     =20.e0 
          tauGABA['SupLTS']to_supaxax     =20.e0 
          tauGABA['SupLTS']to_supLTS      =20.e0 
          tauGABA['SupLTS']to_spinstell   =20.e0 
          tauGABA['SupLTS']to_tuftIB      =20.e0 
          tauGABA['SupLTS']to_tuftRS      =20.e0 
          tauGABA['SupLTS']to_deepbask    =20.e0 
          tauGABA['SupLTS']to_deepaxax    =20.e0 
          tauGABA['SupLTS']to_deepLTS     =20.e0 
          tauGABA['SupLTS']to_nontuftRS   =20.e0  

          tauAMPA_spinstell_to_suppyrRS =2.e0  
          tauNMDA_spinstell_to_suppyrRS =130.e0 
          tauAMPA_spinstell_to_suppyrFRB=2.e0  
          tauNMDA_spinstell_to_suppyrFRB=130.e0 
          tauAMPA_spinstell_to_supbask  =.8e0  
          tauNMDA_spinstell_to_supbask  =100.e0
          tauAMPA_spinstell_to_supaxax  =.8e0  
          tauNMDA_spinstell_to_supaxax  =100.e0
          tauAMPA_spinstell_to_supLTS   =1.e0  
          tauNMDA_spinstell_to_supLTS   =100.e0
          tauAMPA_spinstell_to_spinstell=2.e0  
          tauNMDA_spinstell_to_spinstell=130.e0 
          tauAMPA_spinstell_to_tuftIB   =2.e0  
          tauNMDA_spinstell_to_tuftIB   =130.e0 
          tauAMPA_spinstell_to_tuftRS   =2.e0  
          tauNMDA_spinstell_to_tuftRS   =130.e0
          tauAMPA_spinstell_to_deepbask =.8e0  
          tauNMDA_spinstell_to_deepbask =100.e0
          tauAMPA_spinstell_to_deepaxax =.8e0  
          tauNMDA_spinstell_to_deepaxax =100.e0
          tauAMPA_spinstell_to_deepLTS  =1.e0  
          tauNMDA_spinstell_to_deepLTS  =100.e0
          tauAMPA_spinstell_to_nontuftRS=2.e0  
          tauNMDA_spinstell_to_nontuftRS=130.e0

          tauAMPA_tuftIB_to_suppyrRS    =2.e0 
          tauNMDA_tuftIB_to_suppyrRS    =130.e0
          tauAMPA_tuftIB_to_suppyrFRB   =2.e0 
          tauNMDA_tuftIB_to_suppyrFRB   =130.e0
          tauAMPA_tuftIB_to_supbask     =.8e0  
          tauNMDA_tuftIB_to_supbask     =100.e0 
          tauAMPA_tuftIB_to_supaxax     =.8e0  
          tauNMDA_tuftIB_to_supaxax     =100.e0 
          tauAMPA_tuftIB_to_supLTS      =1.e0  
          tauNMDA_tuftIB_to_supLTS      =100.e0 
          tauAMPA_tuftIB_to_spinstell   =2.e0   
          tauNMDA_tuftIB_to_spinstell   =130.e0 
          tauAMPA_tuftIB_to_tuftIB      =2.e0  
          tauNMDA_tuftIB_to_tuftIB      =130.e0 
          tauAMPA_tuftIB_to_tuftRS      =2.0e0 
          tauNMDA_tuftIB_to_tuftRS      =130.e0 
          tauAMPA_tuftIB_to_deepbask    =.8e0  
          tauNMDA_tuftIB_to_deepbask    =100.e0 
          tauAMPA_tuftIB_to_deepaxax    =.8e0  
          tauNMDA_tuftIB_to_deepaxax    =100.e0 
          tauAMPA_tuftIB_to_deepLTS     =1.e0  
          tauNMDA_tuftIB_to_deepLTS     =100.e0 
          tauAMPA_tuftIB_to_nontuftRS   =2.0e0 
          tauNMDA_tuftIB_to_nontuftRS   =130.e0 

          tauAMPA_tuftRS_to_suppyrRS    =2.e0 
          tauNMDA_tuftRS_to_suppyrRS    =130.e0
          tauAMPA_tuftRS_to_suppyrFRB   =2.e0 
          tauNMDA_tuftRS_to_suppyrFRB   =130.e0
          tauAMPA_tuftRS_to_supbask     =.8e0  
          tauNMDA_tuftRS_to_supbask     =100.e0 
          tauAMPA_tuftRS_to_supaxax     =.8e0  
          tauNMDA_tuftRS_to_supaxax     =100.e0 
          tauAMPA_tuftRS_to_supLTS      =1.e0  
          tauNMDA_tuftRS_to_supLTS      =100.e0 
          tauAMPA_tuftRS_to_spinstell   =2.e0  
          tauNMDA_tuftRS_to_spinstell   =130.e0 
          tauAMPA_tuftRS_to_tuftIB      =2.e0  
          tauNMDA_tuftRS_to_tuftIB      =130.e0 
          tauAMPA_tuftRS_to_tuftRS      =2.e0  
          tauNMDA_tuftRS_to_tuftRS      =130.e0 
          tauAMPA_tuftRS_to_deepbask    =.8e0  
          tauNMDA_tuftRS_to_deepbask    =100.e0 
          tauAMPA_tuftRS_to_deepaxax    =.8e0  
          tauNMDA_tuftRS_to_deepaxax    =100.e0 
          tauAMPA_tuftRS_to_deepLTS     =1.e0   
          tauNMDA_tuftRS_to_deepLTS     =100.e0 
          tauAMPA_tuftRS_to_nontuftRS   =2.e0  
          tauNMDA_tuftRS_to_nontuftRS   =130.e0 

          tauGABA_deepbask_to_spinstell =6.e0  
          tauGABA_deepbask_to_tuftIB    =6.e0  
          tauGABA_deepbask_to_tuftRS    =6.e0  
          tauGABA_deepbask_to_deepbask  =3.e0  
          tauGABA_deepbask_to_deepaxax  =3.e0  
          tauGABA_deepbask_to_deepLTS   =3.e0  
          tauGABA_deepbask_to_nontuftRS =6.e0  

          tauGABA_deepaxax_to_suppyrRS   =6.e0  
          tauGABA_deepaxax_to_suppyrFRB  =6.e0  
          tauGABA_deepaxax_to_spinstell  =6.e0  
          tauGABA_deepaxax_to_tuftIB     =6.e0  
          tauGABA_deepaxax_to_tuftRS     =6.e0  
          tauGABA_deepaxax_to_nontuftRS  =6.e0  

          tauGABA_deepLTS_to_suppyrRS    =20.e0 
          tauGABA_deepLTS_to_suppyrFRB   =20.e0 
          tauGABA_deepLTS_to_supbask     =20.e0 
          tauGABA_deepLTS_to_supaxax     =20.e0 
          tauGABA_deepLTS_to_supLTS      =20.e0 
          tauGABA_deepLTS_to_spinstell   =20.e0 
          tauGABA_deepLTS_to_tuftIB      =20.e0 
          tauGABA_deepLTS_to_tuftRS      =20.e0 
          tauGABA_deepLTS_to_deepbask    =20.e0 
          tauGABA_deepLTS_to_deepaxax    =20.e0 
          tauGABA_deepLTS_to_deepLTS     =20.e0 
          tauGABA_deepLTS_to_nontuftRS   =20.e0 

          tauAMPA_TCR_to_suppyrRS        =2.e0  
          tauNMDA_TCR_to_suppyrRS        =130.e0
          tauAMPA_TCR_to_suppyrFRB       =2.e0  
          tauNMDA_TCR_to_suppyrFRB       =130.e0
          tauAMPA_TCR_to_supbask         =1.e0  
          tauNMDA_TCR_to_supbask         =100.e0
          tauAMPA_TCR_to_supaxax         =1.e0  
          tauNMDA_TCR_to_supaxax         =100.e0 
          tauAMPA_TCR_to_spinstell       =2.0e0 
          tauNMDA_TCR_to_spinstell       =130.e0
          tauAMPA_TCR_to_tuftIB          =2.e0  
          tauNMDA_TCR_to_tuftIB          =130.e0
          tauAMPA_TCR_to_tuftRS          =2.e0  
          tauNMDA_TCR_to_tuftRS          =130.e0
          tauAMPA_TCR_to_deepbask        =1.e0  
          tauNMDA_TCR_to_deepbask        =100.e0
          tauAMPA_TCR_to_deepaxax        =1.e0  
          tauNMDA_TCR_to_deepaxax        =100.e0
          tauAMPA_TCR_to_nRT             =2.0e0      
          tauNMDA_TCR_to_nRT             =150.e0
          tauAMPA_TCR_to_nontuftRS       =2.0e0     
          tauNMDA_TCR_to_nontuftRS       =130.e0

#         tauGABA1_nRT_to_TCR             =10.e0 
#         tauGABA2_nRT_to_TCR             =30.e0 
#         tauGABA1_nRT_to_nRT             =18.e0 
#         tauGABA2_nRT_to_nRT             =89.e0 
# See notebook entry of 17 Feb. 2004.
# Speed these up per Huntsman & Huguenard (2000)
          tauGABA1_nRT_to_TCR             =3.30e0 
          tauGABA2_nRT_to_TCR             =10.e0 
          tauGABA1_nRT_to_nRT             = 9.e0 
          tauGABA2_nRT_to_nRT             =44.5e0 

          tauAMPA_nontuftRS_to_suppyrRS  =2.e0  
          tauNMDA_nontuftRS_to_suppyrRS  =130.e0
          tauAMPA_nontuftRS_to_suppyrFRB =2.e0  
          tauNMDA_nontuftRS_to_suppyrFRB =130.e0
          tauAMPA_nontuftRS_to_supbask   =.8e0  
          tauNMDA_nontuftRS_to_supbask   =100.e0
          tauAMPA_nontuftRS_to_supaxax   =.8e0  
          tauNMDA_nontuftRS_to_supaxax   =100.e0 
          tauAMPA_nontuftRS_to_supLTS    =1.0e0 
          tauNMDA_nontuftRS_to_supLTS    =100.e0
          tauAMPA_nontuftRS_to_spinstell =2.e0  
          tauNMDA_nontuftRS_to_spinstell =130.e0
          tauAMPA_nontuftRS_to_tuftIB    =2.e0  
          tauNMDA_nontuftRS_to_tuftIB    =130.e0
          tauAMPA_nontuftRS_to_tuftRS    =2.e0  
          tauNMDA_nontuftRS_to_tuftRS    =130.e0
          tauAMPA_nontuftRS_to_deepbask  =.8e0  
          tauNMDA_nontuftRS_to_deepbask  =100.e0
          tauAMPA_nontuftRS_to_deepaxax  =.8e0   
          tauNMDA_nontuftRS_to_deepaxax  =100.e0
     
       tauAMPA_nontuftRS_to_deepLTS   =1.e0  
       tauNMDA_nontuftRS_to_deepLTS   =100.e0
       tauAMPA_nontuftRS_to_TCR       =2.e0  
       tauNMDA_nontuftRS_to_TCR       =130.e0 
       tauAMPA_nontuftRS_to_nRT       =2.0e0 
       tauNMDA_nontuftRS_to_nRT       =100.e0 
       tauAMPA_nontuftRS_to_nontuftRS =2.e0  
       tauNMDA_nontuftRS_to_nontuftRS =130.e0 
# End definition of synaptic time constants

# Synaptic conductance scaling factors.
       gAMPA['SupPyrRS']to_suppyrRS =0.25e-3
       gNMDA['SupPyrRS']to_suppyrRS = 0.025e-3
       gAMPA['SupPyrRS']to_suppyrFRB = 0.25e-3
       gNMDA['SupPyrRS']to_suppyrFRB = 0.025e-3
       gAMPA['SupPyrRS']to_supbask  =3.00e-3
       gNMDA['SupPyrRS']to_supbask  =0.15e-3
       gAMPA['SupPyrRS']to_supaxax  =3.0e-3
       gNMDA['SupPyrRS']to_supaxax  =0.15e-3
       gAMPA['SupPyrRS']to_supLTS   =2.0e-3
       gNMDA['SupPyrRS']to_supLTS   =0.15e-3
       gAMPA['SupPyrRS']to_spinstell = 0.10e-3
       gNMDA['SupPyrRS']to_spinstell = 0.01e-3
       gAMPA['SupPyrRS']to_tuftIB   =0.10e-3
       gNMDA['SupPyrRS']to_tuftIB   =0.01e-3
       gAMPA['SupPyrRS']to_tuftRS   =0.10e-3
       gNMDA['SupPyrRS']to_tuftRS   =0.01e-3
       gAMPA['SupPyrRS']to_deepbask =1.00e-3
       gNMDA['SupPyrRS']to_deepbask =0.10e-3
       gAMPA['SupPyrRS']to_deepaxax =1.00e-3
       gNMDA['SupPyrRS']to_deepaxax =0.10e-3
       gAMPA['SupPyrRS']to_deepLTS  =1.00e-3
       gNMDA['SupPyrRS']to_deepLTS  =0.15e-3
       gAMPA['SupPyrRS']to_nontuftRS = 0.50e-3
       gNMDA['SupPyrRS']to_nontuftRS = 0.05e-3

       gAMPA['SupPyrFRB']to_suppyrRS = 0.25e-3
       gNMDA['SupPyrFRB']to_suppyrRS = 0.025e-3
       gAMPA['SupPyrFRB']to_suppyrFRB = 0.25e-3
       gNMDA['SupPyrFRB']to_suppyrFRB = .025e-3
       gAMPA['SupPyrFRB']to_supbask  =3.00e-3
       gNMDA['SupPyrFRB']to_supbask  =0.10e-3
       gAMPA['SupPyrFRB']to_supaxax  =3.0e-3
       gNMDA['SupPyrFRB']to_supaxax  =0.10e-3
       gAMPA['SupPyrFRB']to_supLTS   =2.0e-3
       gNMDA['SupPyrFRB']to_supLTS   =0.10e-3
       gAMPA['SupPyrFRB']to_spinstell = 0.10e-3
       gNMDA['SupPyrFRB']to_spinstell = 0.01e-3
       gAMPA['SupPyrFRB']to_tuftIB   =0.10e-3
       gNMDA['SupPyrFRB']to_tuftIB   =0.01e-3
       gAMPA['SupPyrFRB']to_tuftRS   =0.10e-3
       gNMDA['SupPyrFRB']to_tuftRS   =0.01e-3
       gAMPA['SupPyrFRB']to_deepbask =1.00e-3
       gNMDA['SupPyrFRB']to_deepbask =0.10e-3
       gAMPA['SupPyrFRB']to_deepaxax =1.00e-3
       gNMDA['SupPyrFRB']to_deepaxax =0.10e-3
       gAMPA['SupPyrFRB']to_deepLTS  =1.00e-3
       gNMDA['SupPyrFRB']to_deepLTS  =0.10e-3
       gAMPA['SupPyrFRB']to_nontuftRS = 0.50e-3
       gNMDA['SupPyrFRB']to_nontuftRS = 0.05e-3

       gGABA['SupBasket']to_suppyrRS   =1.2e-3
       gGABA['SupBasket']to_suppyrFRB  =1.2e-3
       gGABA['SupBasket']to_supbask    =0.2e-3
       gGABA['SupBasket']to_supaxax    =0.2e-3
       gGABA['SupBasket']to_supLTS     =0.5e-3
#      gGABA['SupBasket']to_spinstell  =0.7e-3
       gGABA['SupBasket']to_spinstell  =0.1e-3 # if main inhib. to spinstell from deep int.

       gGABA['SupAxoaxonic']to_suppyrRS   =1.2e-3
       gGABA['SupAxoaxonic']to_suppyrFRB  =1.2e-3
#      gGABA['SupAxoaxonic']to_spinstell  =1.0e-3
       gGABA['SupAxoaxonic']to_spinstell  =0.1e-3 # if main inhib. to spinstell from deep int.
       gGABA['SupAxoaxonic']to_tuftIB     =1.0e-3
       gGABA['SupAxoaxonic']to_tuftRS     =1.0e-3
       gGABA['SupAxoaxonic']to_nontuftRS  =1.0e-3


       gGABA['SupLTS']to_suppyrRS    =.01e-3
       gGABA['SupLTS']to_suppyrFRB   =.01e-3
       gGABA['SupLTS']to_supbask     =.01e-3
       gGABA['SupLTS']to_supaxax     =.01e-3
       gGABA['SupLTS']to_supLTS      =.05e-3
       gGABA['SupLTS']to_spinstell   =.01e-3
       gGABA['SupLTS']to_tuftIB      =.02e-3
       gGABA['SupLTS']to_tuftRS      =.02e-3
       gGABA['SupLTS']to_deepbask    =.01e-3
       gGABA['SupLTS']to_deepaxax    =.01e-3
       gGABA['SupLTS']to_deepLTS     =.05e-3
       gGABA['SupLTS']to_nontuftRS   =.01e-3

       gAMPA_spinstell_to_suppyrRS =1.0e-3
       gNMDA_spinstell_to_suppyrRS =0.1e-3
       gAMPA_spinstell_to_suppyrFRB = 1.0e-3
       gNMDA_spinstell_to_suppyrFRB = 0.1e-3
       gAMPA_spinstell_to_supbask  =1.0e-3
       gNMDA_spinstell_to_supbask  =.15e-3
       gAMPA_spinstell_to_supaxax  =1.0e-3
       gNMDA_spinstell_to_supaxax  =.15e-3
       gAMPA_spinstell_to_supLTS   =1.0e-3
       gNMDA_spinstell_to_supLTS   =.15e-3
       gAMPA_spinstell_to_spinstell = 1.0e-3
       gNMDA_spinstell_to_spinstell = 0.1e-3
       gAMPA_spinstell_to_tuftIB   =1.0e-3
       gNMDA_spinstell_to_tuftIB   =0.1e-3
       gAMPA_spinstell_to_tuftRS   =1.0e-3
       gNMDA_spinstell_to_tuftRS   =0.1e-3
       gAMPA_spinstell_to_deepbask =1.0e-3
       gNMDA_spinstell_to_deepbask =.15e-3
       gAMPA_spinstell_to_deepaxax =1.0e-3
       gNMDA_spinstell_to_deepaxax =.15e-3
       gAMPA_spinstell_to_deepLTS  =1.0e-3
       gNMDA_spinstell_to_deepLTS  =.15e-3
       gAMPA_spinstell_to_nontuftRS = 1.0e-3
       gNMDA_spinstell_to_nontuftRS = 0.1e-3

       gAMPA_tuftIB_to_suppyrRS    =0.5e-3
       gNMDA_tuftIB_to_suppyrRS    =0.05e-3
       gAMPA_tuftIB_to_suppyrFRB   =0.5e-3
       gNMDA_tuftIB_to_suppyrFRB   =0.05e-3
       gAMPA_tuftIB_to_supbask     =1.0e-3
       gNMDA_tuftIB_to_supbask     =0.15e-3
       gAMPA_tuftIB_to_supaxax     =1.0e-3
       gNMDA_tuftIB_to_supaxax     =0.15e-3
       gAMPA_tuftIB_to_supLTS      =1.0e-3
       gNMDA_tuftIB_to_supLTS      =0.15e-3
       gAMPA_tuftIB_to_spinstell   =0.5e-3
       gNMDA_tuftIB_to_spinstell   =0.05e-3
       gAMPA_tuftIB_to_tuftIB      =2.0e-3
       gNMDA_tuftIB_to_tuftIB      =0.20e-3
       gAMPA_tuftIB_to_tuftRS      =2.0e-3
       gNMDA_tuftIB_to_tuftRS      =0.20e-3
       gAMPA_tuftIB_to_deepbask    =3.0e-3
       gNMDA_tuftIB_to_deepbask    =0.15e-3
       gAMPA_tuftIB_to_deepaxax    =3.0e-3
       gNMDA_tuftIB_to_deepaxax    =0.15e-3
       gAMPA_tuftIB_to_deepLTS     =2.0e-3
       gNMDA_tuftIB_to_deepLTS     =0.15e-3
       gAMPA_tuftIB_to_nontuftRS   =2.0e-3
       gNMDA_tuftIB_to_nontuftRS   =0.20e-3

       gAMPA_tuftRS_to_suppyrRS    =0.5e-3
       gNMDA_tuftRS_to_suppyrRS    =0.05e-3
       gAMPA_tuftRS_to_suppyrFRB   =0.5e-3
       gNMDA_tuftRS_to_suppyrFRB   =0.05e-3
       gAMPA_tuftRS_to_supbask     =1.0e-3
       gNMDA_tuftRS_to_supbask     =0.15e-3
       gAMPA_tuftRS_to_supaxax     =1.0e-3
       gNMDA_tuftRS_to_supaxax     =0.15e-3
       gAMPA_tuftRS_to_supLTS      =1.0e-3
       gNMDA_tuftRS_to_supLTS      =0.15e-3
       gAMPA_tuftRS_to_spinstell   =0.5e-3
       gNMDA_tuftRS_to_spinstell   =0.05e-3
       gAMPA_tuftRS_to_tuftIB      =1.0e-3
       gNMDA_tuftRS_to_tuftIB      =0.10e-3
       gAMPA_tuftRS_to_tuftRS      =1.0e-3
       gNMDA_tuftRS_to_tuftRS      =0.10e-3
       gAMPA_tuftRS_to_deepbask    =3.0e-3
       gNMDA_tuftRS_to_deepbask    =0.10e-3
       gAMPA_tuftRS_to_deepaxax    =3.0e-3
       gNMDA_tuftRS_to_deepaxax    =0.10e-3
       gAMPA_tuftRS_to_deepLTS     =2.0e-3
       gNMDA_tuftRS_to_deepLTS     =0.10e-3
       gAMPA_tuftRS_to_nontuftRS   =1.0e-3
       gNMDA_tuftRS_to_nontuftRS   =0.10e-3

#      gGABA_deepbask_to_spinstell =1.0e-3
       gGABA_deepbask_to_spinstell =1.5e-3 # ? suppress spiny stellate bursts ?
       gGABA_deepbask_to_tuftIB    =0.7e-3
       gGABA_deepbask_to_tuftRS    =0.7e-3
       gGABA_deepbask_to_deepbask  =0.2e-3
       gGABA_deepbask_to_deepaxax  =0.2e-3
       gGABA_deepbask_to_deepLTS   =0.7e-3
       gGABA_deepbask_to_nontuftRS =0.7e-3

       gGABA_deepaxax_to_suppyrRS   =1.0e-3
       gGABA_deepaxax_to_suppyrFRB  =1.0e-3
#      gGABA_deepaxax_to_spinstell  =1.0e-3
       gGABA_deepaxax_to_spinstell  =1.5e-3 # ? suppress spiny stellate bursts ?
       gGABA_deepaxax_to_tuftIB     =1.0e-3
       gGABA_deepaxax_to_tuftRS     =1.0e-3
       gGABA_deepaxax_to_nontuftRS  =1.0e-3

       gGABA_deepLTS_to_suppyrRS    =.01e-3
       gGABA_deepLTS_to_suppyrFRB   =.01e-3
       gGABA_deepLTS_to_supbask     =.01e-3
       gGABA_deepLTS_to_supaxax     =.01e-3
       gGABA_deepLTS_to_supLTS      =.05e-3
       gGABA_deepLTS_to_spinstell   =.01e-3
#      gGABA_deepLTS_to_tuftIB      =.02e-3
       gGABA_deepLTS_to_tuftIB      =.05e-3 # will this help suppress bursting?
       gGABA_deepLTS_to_tuftRS      =.02e-3
       gGABA_deepLTS_to_deepbask    =.01e-3
       gGABA_deepLTS_to_deepaxax    =.01e-3
       gGABA_deepLTS_to_deepLTS     =.05e-3
       gGABA_deepLTS_to_nontuftRS   =.01e-3

       gAMPA_TCR_to_suppyrRS        =0.5e-3
       gNMDA_TCR_to_suppyrRS        =0.05e-3
       gAMPA_TCR_to_suppyrFRB       =0.5e-3
       gNMDA_TCR_to_suppyrFRB       =0.05e-3
#      gAMPA_TCR_to_supbask         =1.0e-3
       gAMPA_TCR_to_supbask         =0.1e-3
# try a variation in which main feedforward inhibtion from thalamus
# is via deep interneurons.  May be necessary later to include special
# layer 4 interneurons
#      gNMDA_TCR_to_supbask         =.10e-3
       gNMDA_TCR_to_supbask         =.01e-3
#      gAMPA_TCR_to_supaxax         =1.0e-3
       gAMPA_TCR_to_supaxax         =0.1e-3
#      gNMDA_TCR_to_supaxax         =.10e-3
       gNMDA_TCR_to_supaxax         =.01e-3
       gAMPA_TCR_to_spinstell       =1.0e-3
       gNMDA_TCR_to_spinstell       =.10e-3
       gAMPA_TCR_to_tuftIB          =1.5e-3
       gNMDA_TCR_to_tuftIB          =.15e-3
       gAMPA_TCR_to_tuftRS          =1.5e-3
       gNMDA_TCR_to_tuftRS          =.15e-3
#      gAMPA_TCR_to_deepbask        =1.0e-3
       gAMPA_TCR_to_deepbask        =1.5e-3
       gNMDA_TCR_to_deepbask        =.10e-3
       gAMPA_TCR_to_deepaxax        =1.0e-3
       gNMDA_TCR_to_deepaxax        =.10e-3
       gAMPA_TCR_to_nRT             =0.75e-3   
       gNMDA_TCR_to_nRT             =.15e-3
       gAMPA_TCR_to_nontuftRS       =1.0e-3    
       gNMDA_TCR_to_nontuftRS       =.10e-3

#      gGABA_nRT_to_TCR             =1.0e-3
	objref gGABA_nRT_to_TCR
	gGABA_nRT_to_TCR = new Vector(num_nRT+1)

# Values here need to be set below  
       gGABA_nRT_to_nRT             =0.30e-3

       gAMPA_nontuftRS_to_suppyrRS  =0.5e-3
       gNMDA_nontuftRS_to_suppyrRS  =0.05e-3
       gAMPA_nontuftRS_to_suppyrFRB =0.5e-3
       gNMDA_nontuftRS_to_suppyrFRB =0.05e-3
       gAMPA_nontuftRS_to_supbask   =1.0e-3
       gNMDA_nontuftRS_to_supbask   =0.1e-3
       gAMPA_nontuftRS_to_supaxax   =1.0e-3
       gNMDA_nontuftRS_to_supaxax   =0.1e-3
       gAMPA_nontuftRS_to_supLTS    =1.0e-3
       gNMDA_nontuftRS_to_supLTS    =0.1e-3
       gAMPA_nontuftRS_to_spinstell =0.5e-3
       gNMDA_nontuftRS_to_spinstell =0.05e-3
       gAMPA_nontuftRS_to_tuftIB    =1.0e-3
       gNMDA_nontuftRS_to_tuftIB    =0.1e-3
       gAMPA_nontuftRS_to_tuftRS    =1.0e-3
       gNMDA_nontuftRS_to_tuftRS    =0.1e-3
       gAMPA_nontuftRS_to_deepbask  =3.0e-3
       gNMDA_nontuftRS_to_deepbask  =.10e-3
       gAMPA_nontuftRS_to_deepaxax  =3.0e-3
       gNMDA_nontuftRS_to_deepaxax  =.10e-3
       gAMPA_nontuftRS_to_deepLTS   =2.0e-3
       gNMDA_nontuftRS_to_deepLTS   =.10e-3
       gAMPA_nontuftRS_to_TCR       =.75e-3
       gNMDA_nontuftRS_to_TCR       =.075e-3
       gAMPA_nontuftRS_to_nRT       =0.5e-3
       gNMDA_nontuftRS_to_nRT       =0.05e-3
       gAMPA_nontuftRS_to_nontuftRS =1.0e-3
       gNMDA_nontuftRS_to_nontuftRS =0.1e-3
# End defining synaptic conductance scaling factors


def write_file():
    
# 
# synapses.py ends here
