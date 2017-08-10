'''
This module contains a logbinning function adapted by James Clough.
'''

import numpy as np
from numba import jit, int_, float_

#@jit
def log_bin(data, bin_start=1., first_bin_width=1., a=2.,
            datatype='int', drop_zeros=True):
    if drop_zeros:
        data = np.array(data)[data!=0]
    else:
        data = np.array(data)
    max_x = float(np.max(data))
    bin_start = float(bin_start)
    first_bin_width = float(first_bin_width)
    a = float(a)

    max_power = np.ceil(np.log(1+((a-1)*max_x)/first_bin_width)/np.log(a) - 1.)
    widths = first_bin_width*np.power(a,np.arange(max_power+1, dtype='float'))
    bins = np.cumsum(np.concatenate([np.array([bin_start]), widths]))

    indices = np.digitize(data, bins[1:])
    counts = np.bincount(indices)/float(data.size)
    bins = bins[:len(counts)+1]
    widths = widths[:len(counts)]

    bin_indices = np.arange(len(bins)-1)
    bins = np.array(bins)
    if datatype == 'float':
        centres = np.sqrt(np.roll(bins, -1)* bins)[:-1]
    elif datatype == 'int':
        centres = np.empty(len(bin_indices))
        widths = np.diff(np.ceil(bins),1)
        for i in bin_indices:
            centres[i] = geometric_mean(np.arange(np.ceil(bins[i]),
                                        np.ceil(bins[i+1])))
    counts /= widths
    return centres, counts

@jit(float_[:](int_[:]))
def geometric_mean(x):
    s = len(x)
    y = np.log(x)
    z = np.sum(y)/s
    return np.exp(z)
