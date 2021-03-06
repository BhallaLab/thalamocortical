# allowedcomp.py --- 
# 
# Filename: allowedcomp.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Wed Feb 17 22:14:29 2010 (+0530)
# Version: 
# Last-Updated: Mon Jun 28 19:03:06 2010 (+0530)
#           By: Subhasis Ray
#     Update #: 43
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
# This program is free software; you can redistribute it and[or
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

# This horrendous dictionary stores the allowed compartments for
# synaptic connections between pairs of cell type.
ALLOWED_COMP = defaultdict(dict)

ALLOWED_COMP['SupPyrRS']['SupPyrRS'] = [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
					27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
					34,35,36,37]

ALLOWED_COMP['SupPyrRS']['SupPyrFRB'] = [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
					 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
					 34,35,36,37]

ALLOWED_COMP['SupPyrRS']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					 44,45,46,47,48,49]

ALLOWED_COMP['SupPyrRS']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]

ALLOWED_COMP['SupPyrRS']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				      44,45,46,47,48,49]

ALLOWED_COMP['SupPyrRS']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]

ALLOWED_COMP['SupPyrRS']['TuftedIB'] = [39,40,41,42,43,44,45,46]

ALLOWED_COMP['SupPyrRS']['TuftedRS'] = [39,40,41,42,43,44,45,46]

ALLOWED_COMP['SupPyrRS']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]
ALLOWED_COMP['SupPyrRS']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['SupPyrRS']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				       44,45,46,47,48,49]
ALLOWED_COMP['SupPyrRS']['NontuftedRS'] = [38,39,40,41,42,43,44]

ALLOWED_COMP['SupPyrFRB']['SupPyrRS'] = [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
					 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
					 34,35,36,37]
ALLOWED_COMP['SupPyrFRB']['SupPyrFRB'] = [2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
					  27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
					  34,35,36,37]
ALLOWED_COMP['SupPyrFRB']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				       44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					      44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['TuftedIB'] = [39,40,41,42,43,44,45,46]
ALLOWED_COMP['SupPyrFRB']['TuftedRS'] = [39,40,41,42,43,44,45,46]
ALLOWED_COMP['SupPyrFRB']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					   44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					      44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					44,45,46,47,48,49]
ALLOWED_COMP['SupPyrFRB']['NontuftedRS'] = [38,39,40,41,42,43,44]

ALLOWED_COMP['SupBasket']['SupPyrRS'] = [1,2,3,4,5,6,7,8,9,38,39]
ALLOWED_COMP['SupBasket']['SupPyrFRB'] = [1,2,3,4,5,6,7,8,9,38,39]
ALLOWED_COMP['SupBasket']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]
ALLOWED_COMP['SupBasket']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['SupBasket']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				       44,45,46,47,48,49]
ALLOWED_COMP['SupBasket']['SpinyStellate'] = [1,2,15,28,41]

ALLOWED_COMP['SupLTS']['SupPyrRS'] = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
				      31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
				      50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
				      67,68]
ALLOWED_COMP['SupLTS']['SupPyrFRB'] = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
				       31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
				       50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
				       67,68]
ALLOWED_COMP['SupLTS']['SupBasket'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
				       26,27,31,32,33,34,35,36,37,38,39,40,
				       44,45,46,47,48,49,50,51,52,53]
ALLOWED_COMP['SupLTS']['SupAxoaxonic'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
					  26,27,31,32,33,34,35,36,37,38,39,40,
					  44,45,46,47,48,49,50,51,52,53]
ALLOWED_COMP['SupLTS']['SupLTS'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
				    26,27,31,32,33,34,35,36,37,38,39,40,
				    44,45,46,47,48,49,50,51,52,53]
ALLOWED_COMP['SupLTS']['SpinyStellate'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
					   26,27,31,32,33,34,35,36,37,38,39,40,
					   44,45,46,47,48,49,50,51,52,53]
ALLOWED_COMP['SupLTS']['TuftedIB'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
				       29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
				       48,49,50,51,52,53,54,55]
ALLOWED_COMP['SupLTS']['TuftedRS'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
				       29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
				       48,49,50,51,52,53,54,55]
ALLOWED_COMP['SupLTS']['DeepBasket'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
					 47,48,49,50,51]
ALLOWED_COMP['SupLTS']['DeepAxoaxonic'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
					    47,48,49,50,51]
ALLOWED_COMP['SupLTS']['DeepLTS'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
				      47,48,49,50,51]
ALLOWED_COMP['SupLTS']['NontuftedRS'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
					  29,30,31,32,33,34,38,39,40,41,42,43,44]
ALLOWED_COMP['SpinyStellate']['SupPyrRS'] = [ 2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,
					      26,27,28,29,30,31,32,33]
ALLOWED_COMP['SpinyStellate']['SupPyrFRB'] = [ 2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,
					       26,27,28,29,30,31,32,33]
ALLOWED_COMP['SpinyStellate']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					      44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
						 44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					   44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
						  44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['TuftedIB'] = [ 7,8,9,10,11,12,36,37,38,39,40,41]
ALLOWED_COMP['SpinyStellate']['TuftedRS'] = [ 7,8,9,10,11,12,36,37,38,39,40,41]
ALLOWED_COMP['SpinyStellate']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					       44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
						  44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]
ALLOWED_COMP['SpinyStellate']['NontuftedRS'] = [ 37,38,39,40,41]

ALLOWED_COMP['TuftedIB']['SupPyrRS'] = [ 40,41,42,43,44,45,46,47,48,49,50,51,52]
ALLOWED_COMP['TuftedIB']['SupPyrFRB'] = [ 40,41,42,43,44,45,46,47,48,49,50,51,52]
ALLOWED_COMP['TuftedIB']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					 44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				      44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['TuftedIB'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					 38,39,40,41,42,43,44,45,46,47]
ALLOWED_COMP['TuftedIB']['TuftedRS'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					 38,39,40,41,42,43,44,45,46,47]
ALLOWED_COMP['TuftedIB']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				       44,45,46,47,48,49]
ALLOWED_COMP['TuftedIB']['NontuftedRS'] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
					   37,38,39,40,41,42,43,44]

ALLOWED_COMP['TuftedRS']['SupPyrRS'] = [ 40,41,42,43,44,45,46,47,48,49,50,51,52]
ALLOWED_COMP['TuftedRS']['SupPyrFRB'] = [ 40,41,42,43,44,45,46,47,48,49,50,51,52]
ALLOWED_COMP['TuftedRS']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					 44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				      44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['TuftedIB'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					 38,39,40,41,42,43,44,45,46,47]
ALLOWED_COMP['TuftedRS']['TuftedRS'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					 38,39,40,41,42,43,44,45,46,47]
ALLOWED_COMP['TuftedRS']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
				       44,45,46,47,48,49]
ALLOWED_COMP['TuftedRS']['NontuftedRS'] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
					   37,38,39,40,41,42,43,44]

ALLOWED_COMP['DeepBasket']['SpinyStellate'] = [1,2,15,28,41]
ALLOWED_COMP['DeepBasket']['TuftedIB'] = [ 1,2,3,4,5,6,35,36]
ALLOWED_COMP['DeepBasket']['TuftedRS'] = [ 1,2,3,4,5,6,35,36]
ALLOWED_COMP['DeepBasket']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]
ALLOWED_COMP['DeepBasket']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					       44,45,46,47,48,49]
ALLOWED_COMP['DeepBasket']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					 44,45,46,47,48,49]
ALLOWED_COMP['DeepBasket']['NontuftedRS'] = [1,2,3,4,5,6,35,36]

ALLOWED_COMP['DeepLTS']['SupPyrRS'] = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
				       31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
				       50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
				       67,68]
ALLOWED_COMP['DeepLTS']['SupPyrFRB'] = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
					31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
					50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
					67,68]
ALLOWED_COMP['DeepLTS']['SupBasket'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
					 47,48,49,50,51]
ALLOWED_COMP['DeepLTS']['SupAxoaxonic'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
					    47,48,49,50,51] 
ALLOWED_COMP['DeepLTS']['SupLTS'] = [ 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
				      47,48,49,50,51] 
ALLOWED_COMP['DeepLTS']['SpinyStellate'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
					    26,27,31,32,33,34,35,36,37,38,39,40,
					    44,45,46,47,48,49,50,51,52,53]

ALLOWED_COMP['DeepLTS']['TuftedIB'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
					29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
					48,49,50,51,52,53,54,55]

ALLOWED_COMP['DeepLTS']['TuftedRS'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
					29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
					48,49,50,51,52,53,54,55]

ALLOWED_COMP['DeepLTS']['DeepBasket'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
					 26,27,31,32,33,34,35,36,37,38,39,40,
					 44,45,46,47,48,49,50,51,52,53]

ALLOWED_COMP['DeepLTS']['DeepAxoaxonic'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
					    26,27,31,32,33,34,35,36,37,38,39,40,
					    44,45,46,47,48,49,50,51,52,53]

ALLOWED_COMP['DeepLTS']['DeepLTS'] = [5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
				      26,27,31,32,33,34,35,36,37,38,39,40,
				      44,45,46,47,48,49,50,51,52,53]

ALLOWED_COMP['DeepLTS']['NontuftedRS'] = [ 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
					   29,30,31,32,33,34,38,39,40,41,42,43,44] 

ALLOWED_COMP['TCR']['SupPyrRS'] = [45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
				   61,62,63,64,65,66,67,68]

ALLOWED_COMP['TCR']['SupPyrFRB'] = [45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
				    61,62,63,64,65,66,67,68]

ALLOWED_COMP['TCR']['SupBasket'] = [2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['TCR']['SupAxoaxonic'] = [2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['TCR']['SpinyStellate'] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
					37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53]

ALLOWED_COMP['TCR']['TuftedIB'] = [ 47,48,49,50,51,52,53,54,55]

ALLOWED_COMP['TCR']['TuftedRS'] = [ 47,48,49,50,51,52,53,54,55]

ALLOWED_COMP['TCR']['DeepBasket'] = [2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['TCR']['DeepAxoaxonic'] = [2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['TCR']['nRT'] = [2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['TCR']['NontuftedRS'] = [40,41,42,43,44]


ALLOWED_COMP['nRT']['TCR'] = [1,2,15,28,41,54,67,80,93,106,119]

ALLOWED_COMP['nRT']['nRT'] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
			      20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
			      36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,
			      52,53]


ALLOWED_COMP['NontuftedRS']['SupPyrRS'] = [ 41,42,43,44 ]

ALLOWED_COMP['NontuftedRS']['SupPyrFRB'] = [ 41,42,43,44 ]

ALLOWED_COMP['NontuftedRS']['SupBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					    44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['SupAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					       44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['SupLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					 44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['SpinyStellate'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
						44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['TuftedIB'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					    21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					    38,39,40,41,42,43,44,45,46,47]

ALLOWED_COMP['NontuftedRS']['TuftedRS'] = [ 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					    21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
					    38,39,40,41,42,43,44,45,46,47]

ALLOWED_COMP['NontuftedRS']['DeepBasket'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					     44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['DeepAxoaxonic'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
						44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['DeepLTS'] = [5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
					  44,45,46,47,48,49]

ALLOWED_COMP['NontuftedRS']['TCR'] = [ 6, 7, 8, 9, 10, 11, 12, 13, 14,
				       19, 20, 21, 22, 23, 24, 25, 26, 27,
				       32, 33, 34, 35, 36, 37, 38, 39, 40,
				       45, 46, 47, 48, 49, 50, 51, 52, 53,
				       58, 59, 60, 61, 62, 63, 64, 65, 66,
				       71, 72, 73, 74, 75, 76, 77, 78, 79,
				       84, 85, 86, 87, 88, 89, 90, 91, 92,
				       97, 98, 99,100,101,102,103,104,105,
				       110,111,112,113,114,115,116,117,118,
				       123,124,125,126,127,128,129,130,131]

ALLOWED_COMP['NontuftedRS']['nRT'] = [ 2,3,4,15,16,17,28,29,30,41,42,43]

ALLOWED_COMP['NontuftedRS']['NontuftedRS'] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
					      21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
					      37,38,39,40,41,42,43,44]

# import netCDF4 as netcdf
# import numpy
# def write_to_netcdf(filename):
#     dataset = netcdf.Dataset('allowedcompmap.nc', 'w')
#     celltype = numpy.dtype([()])

import csv
if __name__ == '__main__':
    filename = 'allowedcomp.txt'
    outfile = open(filename, 'w')
    output = csv.writer(outfile)
    for pre, entry in ALLOWED_COMP.items():
        for post, comps in entry.items():
            row = [pre, post] + comps
            output.writerow(row)
    outfile.close()

#
# allowedcomp.py ends here
