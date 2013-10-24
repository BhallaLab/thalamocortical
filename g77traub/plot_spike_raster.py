#!/usr/bin/env python

import sys
import numpy as np
import pylab as plt

suppyrRS_base = -1
suppyrFRB_base = 999
supbask_base = 1049
supaxax_base = 1139
supLTS_base = 1229
spinstell_base = 1319
tuftIB_base = 1559
tuftRS_base = 2359
nontuftRS_base = 2559
deepbask_base = 3059
deepaxax_base = 3159
deepLTS_base = 3259
TCR_base = 3359
nRT_base = 3459

colors = {
    'suppyrRS': -1,
    'suppyrFRB': 999,
    'supbask': 1049,
    'supaxax': 1139,
    'supLTS': 1229,
    'spinstell': 1319,
    'tuftIB': 1559,
    'tuftRS': 2359,
    'nontuftRS': 2559,
    'deepbask': 3059,
    'deepaxax': 3159,
    'deepLTS': 3259,
    'TCR': 3359,
    'nRT': 3459,
}
if __name__ == '__main__':
    filename = 'data/out.dat'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    data = np.loadtxt(filename)
    spike_times = {}
    spike_times['suppyrRS'] = data[(data[:,1] > suppyrRS_base) &
                    (data[:,1] <= suppyrFRB_base)]
    spike_times['suppyrFRB'] = data[(data[:,1] > suppyrFRB_base) &
                     (data[:,1] <= supbask_base)]
    spike_times['supbask'] = data[(data[:,1] > supbask_base) &
                   (data[:,1] <= supaxax_base)]
    spike_times['supaxax'] = data[(data[:,1] > supaxax_base) &
                   (data[:,1] <= supLTS_base)]
    spike_times['supLTS'] = data[(data[:,1] > supLTS_base) &
                  (data[:,1] <= spinstell_base)]
    spike_times['spinstell'] = data[(data[:,1] > spinstell_base) &
                     (data[:,1] <= tuftIB_base)]
    spike_times['tuftIB'] = data[(data[:,1] > tuftIB_base) &
                  (data[:,1] <= tuftRS_base)]
    spike_times['tuftRS'] = data[(data[:,1] > tuftRS_base) &
                  (data[:,1] <= nontuftRS_base)]
    spike_times['nontuftRS'] = data[(data[:,1] > nontuftRS_base) &
                     (data[:,1] <= deepbask_base)]
    spike_times['deepbask'] = data[(data[:,1] > deepbask_base) &
                    (data[:,1] <= deepaxax_base)]
    spike_times['deepaxax'] = data[(data[:,1] > deepaxax_base) &
                    (data[:,1] <= deepLTS_base)]
    spike_times['deepLTS'] = data[(data[:,1] > deepLTS_base) &
                   (data[:,1] <= TCR_base)]
    spike_times['TCR'] =  data[(data[:,1] > TCR_base) &
                (data[:,1] <= nRT_base)]
    spike_times['nRT'] =  data[data[:,1] > nRT_base]
    cm = plt.get_cmap('Accent')
    llines = []
    llabels = []
    for ii, (celltype, spikes) in enumerate(spike_times.items()):
        c = cm(float(ii)/len(spike_times))
        llines.append(plt.Line2D([0], [0], color=c, ls='', marker='|', mew=3, mec=c))
        llabels.append(celltype)
        plt.plot(spikes[:, 0], spikes[:, 1], '|', color=c, label=celltype)
    plt.legend(llines, llabels, loc='best')
    plt.show()
