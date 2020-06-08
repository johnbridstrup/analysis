from utilities import getData, half_time, sigmoid, half_index
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

def rangeIdx(initial, final, M):
    dataMax = M[-1]
    idx=0
    while M[idx]/dataMax < initial:
        idx+=1
    idx_i = idx
    while M[idx]/dataMax < final:
        idx+=1
    idx_f = idx
    return idx_i,idx_f

#folder = '/Users/john/Development/Dissertation/analysis/data/k2a_sweep'
folder = '/Users/john/Development/Dissertation/analysis/data/k2a_sweep/expanded/'
data = getData(folder,'k2',runs=True)

plt.figure()
for key,sample in data.items():
    R = math.pow(sample['conc'],sample['n2']-1)*sample['k2']/sample['a']
    # fit = CubicSpline(sample['t'],sample['M'])
    av_ht,av_slope = half_time(sample['M'],sample['t'],True)
    av_run_slope = 0
    sq_run_slope = 0
    for subkey,run in sample['runs'].items():
        ht,slope = half_time(run['M'],run['t'],True)
        # idx1,idx2 = rangeIdx(0.2,0.8,run['M'])
        # del_t = av_ht - ht
        # rsq = 0
        # for idx,t in enumerate(run['t']):
        #     rsq += math.pow(fit(t) - run['M'][idx],2)
        av_run_slope += slope
        sq_run_slope += slope*slope 
    av_run_slope = av_run_slope / sample['ind_runs']
    kok = av_run_slope / av_slope
    plt.scatter(R,kok)
plt.show()
    
