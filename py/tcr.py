# tcr.py --- 
# 
# Filename: tcr.py
# Description: 
# Author: subhasis ray
# Maintainer: 
# Created: Wed May 13 11:28:02 2009 (+0530)
# Version: 
# Last-Updated: Wed Jun 10 00:34:29 2009 (+0530)
#           By: subhasis ray
#     Update #: 44
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

import moose

from compartment import *

class TCR(moose.Cell):
    EK =  -95e-3
    Em =  -70e-3
    ENa =   50e-3
    ECa =   125e-3
    EAR =  -35e-3
    conductance = {
        0: {
            'NaF_TCR':   0.4,
            'KDR':   0.4,
            'KA':   0.001,
            'K2':   0.0005,
            },
        1: {
            'NaF_TCR':   0.1,
            'NaPF_TCR':   0.0002,
            'KDR':   0.075,
            'KC':   0.012,
            'KA':   0.03,
            'KM':   0.0005,
            'K2':   0.002,
            'KAHP_SLOWER':   5.'E'-05,
            'CaL':   0.0005,
            'CaT':   0.0005,
            'AR':   0.00025,
            },
        2: {
            'NaF_TCR':   0.1,
            'NaPF_TCR':   0.0002,
            'KDR':   0.05,
            'KC':   0.012,
            'KA':   0.03,
            'KM':   0.0005,
            'K2':   0.002,
            'KAHP_SLOWER':   5.'E'-05,
            'CaL':   0.0005,
            'CaT':   0.005,
            'AR':   0.0005,
            },
        3: {
            'NaF_TCR':   0.005,
            'NaPF_TCR':   1.'E'-05,
            'KC':   0.02,
            'KA':   0.001,
            'KM':   0.0005,
            'K2':   0.002,
            'KAHP_SLOWER':   5.'E'-05,
            'CaL':   0.00025,
            'CaT':   0.003,
            'AR':   0.0003,
            },
        4: {
            'NaF_TCR':   0.005,
            'NaPF_TCR':   1.'E'-05,
            'KC':   0.02,
            'KA':   0.001,
            'KM':   0.0005,
            'K2':   0.002,
            'KAHP_SLOWER':   5.'E'-05,
            'CaLm':   0.00025,
            'CaT':   0.0005,
            'AR':   0.0003,
            }
        }

       channels = ['AR',
                   'CaL',
                   'CaT',
                   'K2',
                   'KA',
                   'KAHP_SLOWER',
                   'KC',
                   'KDR',
                   'KM',
                   'NaF_TCR',
                   'NaPF_TCR']

       channel_lib = {}
       channel_lib['AR'
       channel_lib['CaL'
       channel_lib['CaT'
       channel_lib['K2'
       channel_lib['KA'
       channel_lib['KAHP_SLOWER'
       channel_lib['KC'
       channel_lib['KDR'
       channel_lib['KM'
       channel_lib['NaF_TCR'
       channel_lib['NaPF_TCR'        
    def __init__(self, *args):
	moose.Cell.__init__(self, *args)
	self.comp = [None]
	self.num_comp = 137
        level = defaultdict(set)
        self.level = level
	for ii in range(1, self.num_comp + 1):
	    self.comp.append(MyCompartment('comp_' + str(ii), self))
	self.comp[1].traubConnect(self.comp[132])
	self.comp[1].traubConnect(self.comp[2])
	self.comp[1].traubConnect(self.comp[15])
	self.comp[1].traubConnect(self.comp[28])
	self.comp[1].traubConnect(self.comp[41])
	self.comp[1].traubConnect(self.comp[54])
	self.comp[1].traubConnect(self.comp[67])
	self.comp[1].traubConnect(self.comp[80])
	self.comp[1].traubConnect(self.comp[93])
	self.comp[1].traubConnect(self.comp[106])
	self.comp[1].traubConnect(self.comp[119])
	self.comp[2].traubConnect(self.comp[3])
	self.comp[2].traubConnect(self.comp[4])
	self.comp[2].traubConnect(self.comp[5])
	self.comp[3].traubConnect(self.comp[6])
	self.comp[3].traubConnect(self.comp[7])
	self.comp[3].traubConnect(self.comp[8])
	self.comp[4].traubConnect(self.comp[9])
	self.comp[4].traubConnect(self.comp[10])
	self.comp[4].traubConnect(self.comp[11])
	self.comp[5].traubConnect(self.comp[12])
	self.comp[5].traubConnect(self.comp[13])
	self.comp[5].traubConnect(self.comp[14])
	self.comp[6].traubConnect(self.comp[7])
	self.comp[6].traubConnect(self.comp[8])
	self.comp[7].traubConnect(self.comp[8])
	self.comp[9].traubConnect(self.comp[10])
	self.comp[9].traubConnect(self.comp[11])
	self.comp[10].traubConnect(self.comp[11])
	self.comp[12].traubConnect(self.comp[13])
	self.comp[12].traubConnect(self.comp[14])
	self.comp[13].traubConnect(self.comp[14])
	self.comp[15].traubConnect(self.comp[16])
	self.comp[15].traubConnect(self.comp[17])
	self.comp[15].traubConnect(self.comp[18])
	self.comp[16].traubConnect(self.comp[19])
	self.comp[16].traubConnect(self.comp[20])
	self.comp[16].traubConnect(self.comp[21])
	self.comp[17].traubConnect(self.comp[22])
	self.comp[17].traubConnect(self.comp[23])
	self.comp[17].traubConnect(self.comp[24])
	self.comp[18].traubConnect(self.comp[25])
	self.comp[18].traubConnect(self.comp[26])
	self.comp[18].traubConnect(self.comp[27])
	self.comp[19].traubConnect(self.comp[20])
	self.comp[19].traubConnect(self.comp[21])
	self.comp[20].traubConnect(self.comp[21])
	self.comp[22].traubConnect(self.comp[23])
	self.comp[22].traubConnect(self.comp[24])
	self.comp[23].traubConnect(self.comp[24])
	self.comp[25].traubConnect(self.comp[26])
	self.comp[25].traubConnect(self.comp[27])
	self.comp[26].traubConnect(self.comp[27])
	self.comp[28].traubConnect(self.comp[29])
	self.comp[28].traubConnect(self.comp[30])
	self.comp[28].traubConnect(self.comp[31])
	self.comp[29].traubConnect(self.comp[32])
	self.comp[29].traubConnect(self.comp[33])
	self.comp[29].traubConnect(self.comp[34])
	self.comp[30].traubConnect(self.comp[35])
	self.comp[30].traubConnect(self.comp[36])
	self.comp[30].traubConnect(self.comp[37])
	self.comp[31].traubConnect(self.comp[38])
	self.comp[31].traubConnect(self.comp[39])
	self.comp[31].traubConnect(self.comp[40])
	self.comp[32].traubConnect(self.comp[33])
	self.comp[32].traubConnect(self.comp[34])
	self.comp[33].traubConnect(self.comp[34])
	self.comp[35].traubConnect(self.comp[36])
	self.comp[35].traubConnect(self.comp[37])
	self.comp[36].traubConnect(self.comp[37])
	self.comp[38].traubConnect(self.comp[39])
	self.comp[38].traubConnect(self.comp[40])
	self.comp[39].traubConnect(self.comp[40])
	self.comp[41].traubConnect(self.comp[42])
	self.comp[41].traubConnect(self.comp[43])
	self.comp[41].traubConnect(self.comp[44])
	self.comp[42].traubConnect(self.comp[45])
	self.comp[42].traubConnect(self.comp[46])
	self.comp[42].traubConnect(self.comp[47])
	self.comp[43].traubConnect(self.comp[48])
	self.comp[43].traubConnect(self.comp[49])
	self.comp[43].traubConnect(self.comp[50])
	self.comp[44].traubConnect(self.comp[51])
	self.comp[44].traubConnect(self.comp[52])
	self.comp[44].traubConnect(self.comp[53])
	self.comp[45].traubConnect(self.comp[46])
	self.comp[45].traubConnect(self.comp[47])
	self.comp[46].traubConnect(self.comp[47])
	self.comp[48].traubConnect(self.comp[49])
	self.comp[48].traubConnect(self.comp[50])
	self.comp[49].traubConnect(self.comp[50])
	self.comp[51].traubConnect(self.comp[52])
	self.comp[51].traubConnect(self.comp[53])
	self.comp[52].traubConnect(self.comp[53])
	self.comp[54].traubConnect(self.comp[55])
	self.comp[54].traubConnect(self.comp[56])
	self.comp[54].traubConnect(self.comp[57])
	self.comp[55].traubConnect(self.comp[58])
	self.comp[55].traubConnect(self.comp[59])
	self.comp[55].traubConnect(self.comp[60])
	self.comp[56].traubConnect(self.comp[61])
	self.comp[56].traubConnect(self.comp[62])
	self.comp[56].traubConnect(self.comp[63])
	self.comp[57].traubConnect(self.comp[64])
	self.comp[57].traubConnect(self.comp[65])
	self.comp[57].traubConnect(self.comp[66])
	self.comp[58].traubConnect(self.comp[59])
	self.comp[58].traubConnect(self.comp[60])
	self.comp[59].traubConnect(self.comp[60])
	self.comp[61].traubConnect(self.comp[62])
	self.comp[61].traubConnect(self.comp[63])
	self.comp[62].traubConnect(self.comp[63])
	self.comp[64].traubConnect(self.comp[65])
	self.comp[64].traubConnect(self.comp[66])
	self.comp[65].traubConnect(self.comp[66])
	self.comp[67].traubConnect(self.comp[68])
	self.comp[67].traubConnect(self.comp[69])
	self.comp[67].traubConnect(self.comp[70])
	self.comp[68].traubConnect(self.comp[71])
	self.comp[68].traubConnect(self.comp[72])
	self.comp[68].traubConnect(self.comp[73])
	self.comp[69].traubConnect(self.comp[74])
	self.comp[69].traubConnect(self.comp[75])
	self.comp[69].traubConnect(self.comp[76])
	self.comp[70].traubConnect(self.comp[77])
	self.comp[70].traubConnect(self.comp[78])
	self.comp[70].traubConnect(self.comp[79])
	self.comp[71].traubConnect(self.comp[72])
	self.comp[71].traubConnect(self.comp[73])
	self.comp[72].traubConnect(self.comp[73])
	self.comp[74].traubConnect(self.comp[75])
	self.comp[74].traubConnect(self.comp[76])
	self.comp[75].traubConnect(self.comp[76])
	self.comp[77].traubConnect(self.comp[78])
	self.comp[77].traubConnect(self.comp[79])
	self.comp[78].traubConnect(self.comp[79])
	self.comp[80].traubConnect(self.comp[81])
	self.comp[80].traubConnect(self.comp[82])
	self.comp[80].traubConnect(self.comp[83])
	self.comp[81].traubConnect(self.comp[84])
	self.comp[81].traubConnect(self.comp[85])
	self.comp[81].traubConnect(self.comp[86])
	self.comp[82].traubConnect(self.comp[87])
	self.comp[82].traubConnect(self.comp[88])
	self.comp[82].traubConnect(self.comp[89])
	self.comp[83].traubConnect(self.comp[90])
	self.comp[83].traubConnect(self.comp[91])
	self.comp[83].traubConnect(self.comp[92])
	self.comp[84].traubConnect(self.comp[85])
	self.comp[84].traubConnect(self.comp[86])
	self.comp[85].traubConnect(self.comp[86])
	self.comp[87].traubConnect(self.comp[88])
	self.comp[87].traubConnect(self.comp[89])
	self.comp[88].traubConnect(self.comp[89])
	self.comp[90].traubConnect(self.comp[91])
	self.comp[90].traubConnect(self.comp[92])
	self.comp[91].traubConnect(self.comp[92])
	self.comp[93].traubConnect(self.comp[94])
	self.comp[93].traubConnect(self.comp[95])
	self.comp[93].traubConnect(self.comp[96])
	self.comp[94].traubConnect(self.comp[97])
	self.comp[94].traubConnect(self.comp[98])
	self.comp[94].traubConnect(self.comp[99])
	self.comp[95].traubConnect(self.comp[100])
	self.comp[95].traubConnect(self.comp[101])
	self.comp[95].traubConnect(self.comp[102])
	self.comp[96].traubConnect(self.comp[103])
	self.comp[96].traubConnect(self.comp[104])
	self.comp[96].traubConnect(self.comp[105])
	self.comp[97].traubConnect(self.comp[98])
	self.comp[97].traubConnect(self.comp[99])
	self.comp[98].traubConnect(self.comp[99])
	self.comp[100].traubConnect(self.comp[101])
	self.comp[100].traubConnect(self.comp[102])
	self.comp[101].traubConnect(self.comp[102])
	self.comp[103].traubConnect(self.comp[104])
	self.comp[103].traubConnect(self.comp[105])
	self.comp[104].traubConnect(self.comp[105])
	self.comp[106].traubConnect(self.comp[107])
	self.comp[106].traubConnect(self.comp[108])
	self.comp[106].traubConnect(self.comp[109])
	self.comp[107].traubConnect(self.comp[110])
	self.comp[107].traubConnect(self.comp[111])
	self.comp[107].traubConnect(self.comp[112])
	self.comp[108].traubConnect(self.comp[113])
	self.comp[108].traubConnect(self.comp[114])
	self.comp[108].traubConnect(self.comp[115])
	self.comp[109].traubConnect(self.comp[116])
	self.comp[109].traubConnect(self.comp[117])
	self.comp[109].traubConnect(self.comp[118])
	self.comp[110].traubConnect(self.comp[111])
	self.comp[110].traubConnect(self.comp[112])
	self.comp[111].traubConnect(self.comp[112])
	self.comp[113].traubConnect(self.comp[114])
	self.comp[113].traubConnect(self.comp[115])
	self.comp[114].traubConnect(self.comp[115])
	self.comp[116].traubConnect(self.comp[117])
	self.comp[116].traubConnect(self.comp[118])
	self.comp[117].traubConnect(self.comp[118])
	self.comp[119].traubConnect(self.comp[120])
	self.comp[119].traubConnect(self.comp[121])
	self.comp[119].traubConnect(self.comp[122])
	self.comp[120].traubConnect(self.comp[123])
	self.comp[120].traubConnect(self.comp[124])
	self.comp[120].traubConnect(self.comp[125])
	self.comp[121].traubConnect(self.comp[126])
	self.comp[121].traubConnect(self.comp[127])
	self.comp[121].traubConnect(self.comp[128])
	self.comp[122].traubConnect(self.comp[129])
	self.comp[122].traubConnect(self.comp[130])
	self.comp[122].traubConnect(self.comp[131])
	self.comp[123].traubConnect(self.comp[124])
	self.comp[123].traubConnect(self.comp[125])
	self.comp[124].traubConnect(self.comp[125])
	self.comp[126].traubConnect(self.comp[127])
	self.comp[126].traubConnect(self.comp[128])
	self.comp[127].traubConnect(self.comp[128])
	self.comp[129].traubConnect(self.comp[130])
	self.comp[129].traubConnect(self.comp[131])
	self.comp[130].traubConnect(self.comp[131])
	self.comp[132].traubConnect(self.comp[133])
	self.comp[133].traubConnect(self.comp[134])
	self.comp[133].traubConnect(self.comp[136])
	self.comp[134].traubConnect(self.comp[135])
	self.comp[134].traubConnect(self.comp[136])
	self.comp[136].traubConnect(self.comp[137])


 	level[ 1].add(comp[ 1])
 	level[ 2].add(comp[ 2])
 	level[ 3].add(comp[ 3])
 	level[ 3].add(comp[ 4])
 	level[ 3].add(comp[ 5])
 	level[ 4].add(comp[ 6])
 	level[ 4].add(comp[ 7])
 	level[ 4].add(comp[ 8])
 	level[ 4].add(comp[ 9])
 	level[ 4].add(comp[ 10])
 	level[ 4].add(comp[ 11])
 	level[ 4].add(comp[ 12])
 	level[ 4].add(comp[ 13])
 	level[ 4].add(comp[ 14])
 	level[ 2].add(comp[ 15])
 	level[ 3].add(comp[ 16])
 	level[ 3].add(comp[ 17])
 	level[ 3].add(comp[ 18])
 	level[ 4].add(comp[ 19])
 	level[ 4].add(comp[ 20])
 	level[ 4].add(comp[ 21])
 	level[ 4].add(comp[ 22])
 	level[ 4].add(comp[ 23])
 	level[ 4].add(comp[ 24])
 	level[ 4].add(comp[ 25])
 	level[ 4].add(comp[ 26])
 	level[ 4].add(comp[ 27])
 	level[ 2].add(comp[ 28])
 	level[ 3].add(comp[ 29])
 	level[ 3].add(comp[ 30])
 	level[ 3].add(comp[ 31])
 	level[ 4].add(comp[ 32])
 	level[ 4].add(comp[ 33])
 	level[ 4].add(comp[ 34])
 	level[ 4].add(comp[ 35])
 	level[ 4].add(comp[ 36])
 	level[ 4].add(comp[ 37])
 	level[ 4].add(comp[ 38])
 	level[ 4].add(comp[ 39])
 	level[ 4].add(comp[ 40])
 	level[ 2].add(comp[ 41])
 	level[ 3].add(comp[ 42])
 	level[ 3].add(comp[ 43])
 	level[ 3].add(comp[ 44])
 	level[ 4].add(comp[ 45])
 	level[ 4].add(comp[ 46])
 	level[ 4].add(comp[ 47])
 	level[ 4].add(comp[ 48])
 	level[ 4].add(comp[ 49])
 	level[ 4].add(comp[ 50])
 	level[ 4].add(comp[ 51])
 	level[ 4].add(comp[ 52])
 	level[ 4].add(comp[ 53])
 	level[ 2].add(comp[ 54])
 	level[ 3].add(comp[ 55])
 	level[ 3].add(comp[ 56])
 	level[ 3].add(comp[ 57])
 	level[ 4].add(comp[ 58])
 	level[ 4].add(comp[ 59])
 	level[ 4].add(comp[ 60])
 	level[ 4].add(comp[ 61])
 	level[ 4].add(comp[ 62])
 	level[ 4].add(comp[ 63])
 	level[ 4].add(comp[ 64])
 	level[ 4].add(comp[ 65])
 	level[ 4].add(comp[ 66])
 	level[ 2].add(comp[ 67])
 	level[ 3].add(comp[ 68])
 	level[ 3].add(comp[ 69])
 	level[ 3].add(comp[ 70])
 	level[ 4].add(comp[ 71])
 	level[ 4].add(comp[ 72])
 	level[ 4].add(comp[ 73])
 	level[ 4].add(comp[ 74])
 	level[ 4].add(comp[ 75])
 	level[ 4].add(comp[ 76])
 	level[ 4].add(comp[ 77])
 	level[ 4].add(comp[ 78])
 	level[ 4].add(comp[ 79])
 	level[ 2].add(comp[ 80])
 	level[ 3].add(comp[ 81])
 	level[ 3].add(comp[ 82])
 	level[ 3].add(comp[ 83])
 	level[ 4].add(comp[ 84])
 	level[ 4].add(comp[ 85])
 	level[ 4].add(comp[ 86])
 	level[ 4].add(comp[ 87])
 	level[ 4].add(comp[ 88])
 	level[ 4].add(comp[ 89])
 	level[ 4].add(comp[ 90])
 	level[ 4].add(comp[ 91])
 	level[ 4].add(comp[ 92])
 	level[ 2].add(comp[ 93])
 	level[ 3].add(comp[ 94])
 	level[ 3].add(comp[ 95])
 	level[ 3].add(comp[ 96])
 	level[ 4].add(comp[ 97])
 	level[ 4].add(comp[ 98])
 	level[ 4].add(comp[ 99])
 	level[ 4].add(comp[ 100])
 	level[ 4].add(comp[ 101])
 	level[ 4].add(comp[ 102])
 	level[ 4].add(comp[ 103])
 	level[ 4].add(comp[ 104])
 	level[ 4].add(comp[ 105])
 	level[ 2].add(comp[ 106])
 	level[ 3].add(comp[ 107])
 	level[ 3].add(comp[ 108])
 	level[ 3].add(comp[ 109])
 	level[ 4].add(comp[ 110])
 	level[ 4].add(comp[ 111])
 	level[ 4].add(comp[ 112])
 	level[ 4].add(comp[ 113])
 	level[ 4].add(comp[ 114])
 	level[ 4].add(comp[ 115])
 	level[ 4].add(comp[ 116])
 	level[ 4].add(comp[ 117])
 	level[ 4].add(comp[ 118])
 	level[ 2].add(comp[ 119])
 	level[ 3].add(comp[ 120])
 	level[ 3].add(comp[ 121])
 	level[ 3].add(comp[ 122])
 	level[ 4].add(comp[ 123])
 	level[ 4].add(comp[ 124])
 	level[ 4].add(comp[ 125])
 	level[ 4].add(comp[ 126])
 	level[ 4].add(comp[ 127])
 	level[ 4].add(comp[ 128])
 	level[ 4].add(comp[ 129])
 	level[ 4].add(comp[ 130])
 	level[ 4].add(comp[ 131])
 	level[ 0].add(comp[ 132])
 	level[ 0].add(comp[ 133])
 	level[ 0].add(comp[ 134])
 	level[ 0].add(comp[ 135])
 	level[ 0].add(comp[ 136])
 	level[ 0].add(comp[ 137])
	comp[ 1].diameter = 2*  10. * 1e-6
	comp[ 2].diameter = 2*  0.73 * 1e-6
	comp[ 3].diameter = 2*  0.584 * 1e-6
	comp[ 4].diameter = 2*  0.584 * 1e-6
	comp[ 5].diameter = 2*  0.584 * 1e-6
	comp[ 6].diameter = 2*  0.438 * 1e-6
	comp[ 7].diameter = 2*  0.438 * 1e-6
	comp[ 8].diameter = 2*  0.438 * 1e-6
	comp[ 9].diameter = 2*  0.438 * 1e-6
	comp[ 10].diameter = 2*  0.438 * 1e-6
	comp[ 11].diameter = 2*  0.438 * 1e-6
	comp[ 12].diameter = 2*  0.438 * 1e-6
	comp[ 13].diameter = 2*  0.438 * 1e-6
	comp[ 14].diameter = 2*  0.438 * 1e-6
	comp[ 15].diameter = 2*  0.73 * 1e-6
	comp[ 16].diameter = 2*  0.584 * 1e-6
	comp[ 17].diameter = 2*  0.584 * 1e-6
	comp[ 18].diameter = 2*  0.584 * 1e-6
	comp[ 19].diameter = 2*  0.438 * 1e-6
	comp[ 20].diameter = 2*  0.438 * 1e-6
	comp[ 21].diameter = 2*  0.438 * 1e-6
	comp[ 22].diameter = 2*  0.438 * 1e-6
	comp[ 23].diameter = 2*  0.438 * 1e-6
	comp[ 24].diameter = 2*  0.438 * 1e-6
	comp[ 25].diameter = 2*  0.438 * 1e-6
	comp[ 26].diameter = 2*  0.438 * 1e-6
	comp[ 27].diameter = 2*  0.438 * 1e-6
	comp[ 28].diameter = 2*  0.73 * 1e-6
	comp[ 29].diameter = 2*  0.584 * 1e-6
	comp[ 30].diameter = 2*  0.584 * 1e-6
	comp[ 31].diameter = 2*  0.584 * 1e-6
	comp[ 32].diameter = 2*  0.438 * 1e-6
	comp[ 33].diameter = 2*  0.438 * 1e-6
	comp[ 34].diameter = 2*  0.438 * 1e-6
	comp[ 35].diameter = 2*  0.438 * 1e-6
	comp[ 36].diameter = 2*  0.438 * 1e-6
	comp[ 37].diameter = 2*  0.438 * 1e-6
	comp[ 38].diameter = 2*  0.438 * 1e-6
	comp[ 39].diameter = 2*  0.438 * 1e-6
	comp[ 40].diameter = 2*  0.438 * 1e-6
	comp[ 41].diameter = 2*  0.73 * 1e-6
	comp[ 42].diameter = 2*  0.584 * 1e-6
	comp[ 43].diameter = 2*  0.584 * 1e-6
	comp[ 44].diameter = 2*  0.584 * 1e-6
	comp[ 45].diameter = 2*  0.438 * 1e-6
	comp[ 46].diameter = 2*  0.438 * 1e-6
	comp[ 47].diameter = 2*  0.438 * 1e-6
	comp[ 48].diameter = 2*  0.438 * 1e-6
	comp[ 49].diameter = 2*  0.438 * 1e-6
	comp[ 50].diameter = 2*  0.438 * 1e-6
	comp[ 51].diameter = 2*  0.438 * 1e-6
	comp[ 52].diameter = 2*  0.438 * 1e-6
	comp[ 53].diameter = 2*  0.438 * 1e-6
	comp[ 54].diameter = 2*  0.73 * 1e-6
	comp[ 55].diameter = 2*  0.584 * 1e-6
	comp[ 56].diameter = 2*  0.584 * 1e-6
	comp[ 57].diameter = 2*  0.584 * 1e-6
	comp[ 58].diameter = 2*  0.438 * 1e-6
	comp[ 59].diameter = 2*  0.438 * 1e-6
	comp[ 60].diameter = 2*  0.438 * 1e-6
	comp[ 61].diameter = 2*  0.438 * 1e-6
	comp[ 62].diameter = 2*  0.438 * 1e-6
	comp[ 63].diameter = 2*  0.438 * 1e-6
	comp[ 64].diameter = 2*  0.438 * 1e-6
	comp[ 65].diameter = 2*  0.438 * 1e-6
	comp[ 66].diameter = 2*  0.438 * 1e-6
	comp[ 67].diameter = 2*  0.73 * 1e-6
	comp[ 68].diameter = 2*  0.584 * 1e-6
	comp[ 69].diameter = 2*  0.584 * 1e-6
	comp[ 70].diameter = 2*  0.584 * 1e-6
	comp[ 71].diameter = 2*  0.438 * 1e-6
	comp[ 72].diameter = 2*  0.438 * 1e-6
	comp[ 73].diameter = 2*  0.438 * 1e-6
	comp[ 74].diameter = 2*  0.438 * 1e-6
	comp[ 75].diameter = 2*  0.438 * 1e-6
	comp[ 76].diameter = 2*  0.438 * 1e-6
	comp[ 77].diameter = 2*  0.438 * 1e-6
	comp[ 78].diameter = 2*  0.438 * 1e-6
	comp[ 79].diameter = 2*  0.438 * 1e-6
	comp[ 80].diameter = 2*  0.73 * 1e-6
	comp[ 81].diameter = 2*  0.584 * 1e-6
	comp[ 82].diameter = 2*  0.584 * 1e-6
	comp[ 83].diameter = 2*  0.584 * 1e-6
	comp[ 84].diameter = 2*  0.438 * 1e-6
	comp[ 85].diameter = 2*  0.438 * 1e-6
	comp[ 86].diameter = 2*  0.438 * 1e-6
	comp[ 87].diameter = 2*  0.438 * 1e-6
	comp[ 88].diameter = 2*  0.438 * 1e-6
	comp[ 89].diameter = 2*  0.438 * 1e-6
	comp[ 90].diameter = 2*  0.438 * 1e-6
	comp[ 91].diameter = 2*  0.438 * 1e-6
	comp[ 92].diameter = 2*  0.438 * 1e-6
	comp[ 93].diameter = 2*  0.73 * 1e-6
	comp[ 94].diameter = 2*  0.584 * 1e-6
	comp[ 95].diameter = 2*  0.584 * 1e-6
	comp[ 96].diameter = 2*  0.584 * 1e-6
	comp[ 97].diameter = 2*  0.438 * 1e-6
	comp[ 98].diameter = 2*  0.438 * 1e-6
	comp[ 99].diameter = 2*  0.438 * 1e-6
	comp[ 100].diameter = 2*  0.438 * 1e-6
	comp[ 101].diameter = 2*  0.438 * 1e-6
	comp[ 102].diameter = 2*  0.438 * 1e-6
	comp[ 103].diameter = 2*  0.438 * 1e-6
	comp[ 104].diameter = 2*  0.438 * 1e-6
	comp[ 105].diameter = 2*  0.438 * 1e-6
	comp[ 106].diameter = 2*  0.73 * 1e-6
	comp[ 107].diameter = 2*  0.584 * 1e-6
	comp[ 108].diameter = 2*  0.584 * 1e-6
	comp[ 109].diameter = 2*  0.584 * 1e-6
	comp[ 110].diameter = 2*  0.438 * 1e-6
	comp[ 111].diameter = 2*  0.438 * 1e-6
	comp[ 112].diameter = 2*  0.438 * 1e-6
	comp[ 113].diameter = 2*  0.438 * 1e-6
	comp[ 114].diameter = 2*  0.438 * 1e-6
	comp[ 115].diameter = 2*  0.438 * 1e-6
	comp[ 116].diameter = 2*  0.438 * 1e-6
	comp[ 117].diameter = 2*  0.438 * 1e-6
	comp[ 118].diameter = 2*  0.438 * 1e-6
	comp[ 119].diameter = 2*  0.73 * 1e-6
	comp[ 120].diameter = 2*  0.584 * 1e-6
	comp[ 121].diameter = 2*  0.584 * 1e-6
	comp[ 122].diameter = 2*  0.584 * 1e-6
	comp[ 123].diameter = 2*  0.438 * 1e-6
	comp[ 124].diameter = 2*  0.438 * 1e-6
	comp[ 125].diameter = 2*  0.438 * 1e-6
	comp[ 126].diameter = 2*  0.438 * 1e-6
	comp[ 127].diameter = 2*  0.438 * 1e-6
	comp[ 128].diameter = 2*  0.438 * 1e-6
	comp[ 129].diameter = 2*  0.438 * 1e-6
	comp[ 130].diameter = 2*  0.438 * 1e-6
	comp[ 131].diameter = 2*  0.438 * 1e-6
	comp[ 132].diameter = 2*  0.8 * 1e-6
	comp[ 133].diameter = 2*  0.7 * 1e-6
	comp[ 134].diameter = 2*  0.5 * 1e-6
	comp[ 135].diameter = 2*  0.5 * 1e-6
	comp[ 136].diameter = 2*  0.5 * 1e-6
	comp[ 137].diameter = 2*  0.5 * 1e-6
	comp[ 1].length = 42. * 1e-6
	comp[ 2].length = 20. * 1e-6
	comp[ 3].length = 57.5 * 1e-6
	comp[ 4].length = 57.5 * 1e-6
	comp[ 5].length = 57.5 * 1e-6
	comp[ 6].length = 57.5 * 1e-6
	comp[ 7].length = 57.5 * 1e-6
	comp[ 8].length = 57.5 * 1e-6
	comp[ 9].length = 57.5 * 1e-6
	comp[ 10].length = 57.5 * 1e-6
	comp[ 11].length = 57.5 * 1e-6
	comp[ 12].length = 57.5 * 1e-6
	comp[ 13].length = 57.5 * 1e-6
	comp[ 14].length = 57.5 * 1e-6
	comp[ 15].length = 20. * 1e-6
	comp[ 16].length = 57.5 * 1e-6
	comp[ 17].length = 57.5 * 1e-6
	comp[ 18].length = 57.5 * 1e-6
	comp[ 19].length = 57.5 * 1e-6
	comp[ 20].length = 57.5 * 1e-6
	comp[ 21].length = 57.5 * 1e-6
	comp[ 22].length = 57.5 * 1e-6
	comp[ 23].length = 57.5 * 1e-6
	comp[ 24].length = 57.5 * 1e-6
	comp[ 25].length = 57.5 * 1e-6
	comp[ 26].length = 57.5 * 1e-6
	comp[ 27].length = 57.5 * 1e-6
	comp[ 28].length = 20. * 1e-6
	comp[ 29].length = 57.5 * 1e-6
	comp[ 30].length = 57.5 * 1e-6
	comp[ 31].length = 57.5 * 1e-6
	comp[ 32].length = 57.5 * 1e-6
	comp[ 33].length = 57.5 * 1e-6
	comp[ 34].length = 57.5 * 1e-6
	comp[ 35].length = 57.5 * 1e-6
	comp[ 36].length = 57.5 * 1e-6
	comp[ 37].length = 57.5 * 1e-6
	comp[ 38].length = 57.5 * 1e-6
	comp[ 39].length = 57.5 * 1e-6
	comp[ 40].length = 57.5 * 1e-6
	comp[ 41].length = 20. * 1e-6
	comp[ 42].length = 57.5 * 1e-6
	comp[ 43].length = 57.5 * 1e-6
	comp[ 44].length = 57.5 * 1e-6
	comp[ 45].length = 57.5 * 1e-6
	comp[ 46].length = 57.5 * 1e-6
	comp[ 47].length = 57.5 * 1e-6
	comp[ 48].length = 57.5 * 1e-6
	comp[ 49].length = 57.5 * 1e-6
	comp[ 50].length = 57.5 * 1e-6
	comp[ 51].length = 57.5 * 1e-6
	comp[ 52].length = 57.5 * 1e-6
	comp[ 53].length = 57.5 * 1e-6
	comp[ 54].length = 20. * 1e-6
	comp[ 55].length = 57.5 * 1e-6
	comp[ 56].length = 57.5 * 1e-6
	comp[ 57].length = 57.5 * 1e-6
	comp[ 58].length = 57.5 * 1e-6
	comp[ 59].length = 57.5 * 1e-6
	comp[ 60].length = 57.5 * 1e-6
	comp[ 61].length = 57.5 * 1e-6
	comp[ 62].length = 57.5 * 1e-6
	comp[ 63].length = 57.5 * 1e-6
	comp[ 64].length = 57.5 * 1e-6
	comp[ 65].length = 57.5 * 1e-6
	comp[ 66].length = 57.5 * 1e-6
	comp[ 67].length = 20. * 1e-6
	comp[ 68].length = 57.5 * 1e-6
	comp[ 69].length = 57.5 * 1e-6
	comp[ 70].length = 57.5 * 1e-6
	comp[ 71].length = 57.5 * 1e-6
	comp[ 72].length = 57.5 * 1e-6
	comp[ 73].length = 57.5 * 1e-6
	comp[ 74].length = 57.5 * 1e-6
	comp[ 75].length = 57.5 * 1e-6
	comp[ 76].length = 57.5 * 1e-6
	comp[ 77].length = 57.5 * 1e-6
	comp[ 78].length = 57.5 * 1e-6
	comp[ 79].length = 57.5 * 1e-6
	comp[ 80].length = 20. * 1e-6
	comp[ 81].length = 57.5 * 1e-6
	comp[ 82].length = 57.5 * 1e-6
	comp[ 83].length = 57.5 * 1e-6
	comp[ 84].length = 57.5 * 1e-6
	comp[ 85].length = 57.5 * 1e-6
	comp[ 86].length = 57.5 * 1e-6
	comp[ 87].length = 57.5 * 1e-6
	comp[ 88].length = 57.5 * 1e-6
	comp[ 89].length = 57.5 * 1e-6
	comp[ 90].length = 57.5 * 1e-6
	comp[ 91].length = 57.5 * 1e-6
	comp[ 92].length = 57.5 * 1e-6
	comp[ 93].length = 20. * 1e-6
	comp[ 94].length = 57.5 * 1e-6
	comp[ 95].length = 57.5 * 1e-6
	comp[ 96].length = 57.5 * 1e-6
	comp[ 97].length = 57.5 * 1e-6
	comp[ 98].length = 57.5 * 1e-6
	comp[ 99].length = 57.5 * 1e-6
	comp[ 100].length = 57.5 * 1e-6
	comp[ 101].length = 57.5 * 1e-6
	comp[ 102].length = 57.5 * 1e-6
	comp[ 103].length = 57.5 * 1e-6
	comp[ 104].length = 57.5 * 1e-6
	comp[ 105].length = 57.5 * 1e-6
	comp[ 106].length = 20. * 1e-6
	comp[ 107].length = 57.5 * 1e-6
	comp[ 108].length = 57.5 * 1e-6
	comp[ 109].length = 57.5 * 1e-6
	comp[ 110].length = 57.5 * 1e-6
	comp[ 111].length = 57.5 * 1e-6
	comp[ 112].length = 57.5 * 1e-6
	comp[ 113].length = 57.5 * 1e-6
	comp[ 114].length = 57.5 * 1e-6
	comp[ 115].length = 57.5 * 1e-6
	comp[ 116].length = 57.5 * 1e-6
	comp[ 117].length = 57.5 * 1e-6
	comp[ 118].length = 57.5 * 1e-6
	comp[ 119].length = 20. * 1e-6
	comp[ 120].length = 57.5 * 1e-6
	comp[ 121].length = 57.5 * 1e-6
	comp[ 122].length = 57.5 * 1e-6
	comp[ 123].length = 57.5 * 1e-6
	comp[ 124].length = 57.5 * 1e-6
	comp[ 125].length = 57.5 * 1e-6
	comp[ 126].length = 57.5 * 1e-6
	comp[ 127].length = 57.5 * 1e-6
	comp[ 128].length = 57.5 * 1e-6
	comp[ 129].length = 57.5 * 1e-6
	comp[ 130].length = 57.5 * 1e-6
	comp[ 131].length = 57.5 * 1e-6
	comp[ 132].length = 50. * 1e-6
	comp[ 133].length = 50. * 1e-6
	comp[ 134].length = 50. * 1e-6
	comp[ 135].length = 50. * 1e-6
	comp[ 136].length = 50. * 1e-6
	comp[ 137].length = 50. * 1e-6


# 
# tcr.py ends here
