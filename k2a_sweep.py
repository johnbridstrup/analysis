from utilities import getData, half_time, sigmoid, half_index
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit

#folder = '/Users/john/Development/Dissertation/analysis/data/k2a_sweep'
folder = '/Users/john/Development/Dissertation/analysis/data/k2a_sweep/n3'
data = getData(folder,'k2',runs=True)


for key,sample in data.items():
    data[key]['t_half'] = half_time(sample['M'], sample['t'])
    R = math.pow(sample['conc'],sample['n2']-1)*sample['k2']/sample['a']
    plt.figure(1)
    plt.scatter(R, data[key]['t_half'])
    mMax = sample['M'][-1]
    half_guess = sample['t'][half_index(sample['M'])]
    k_guess = mMax/(40*half_guess)
    guess = [mMax, k_guess, half_guess]
    fit,_ = curve_fit(sigmoid, sample['t'], sample['M'], guess)
    data[key]['k'] = fit[1]
    thav=0
    thsq=0
    kav=0
    ksq=0
    kcb=0
    for subkey,run in sample['runs'].items():
        th = half_time(run['M'],run['t'])
        data[key]['runs'][subkey]['t_half'] = th
        thav += th
        thsq += th*th
        mMax = run['M'][-1]
        half_guess = run['t'][half_index(run['M'])]
        k_guess = mMax/(40*half_guess)
        guess = [mMax, k_guess, half_guess]
        fit,_ = curve_fit(sigmoid, run['t'], run['M'], guess)
        kav+=fit[1]
        ksq+=fit[1]*fit[1]
        kcb+=fit[1]*fit[1]*fit[1]
    thav = thav / sample['ind_runs']
    thsq = thsq / sample['ind_runs']
    kav = kav / sample['ind_runs']
    ksq = ksq / sample['ind_runs']
    kcb = kcb / sample['ind_runs']
    k_rel = kav / data[key]['k']
    th_dev = np.sqrt(thsq - (thav*thav))
    th_rel_dev = th_dev / data[key]['t_half']
    data[key]['t_half_dev'] = th_dev
    data[key]['t_half_rel_dev'] = th_rel_dev
    plt.figure(2)
    plt.scatter(R,th_rel_dev)
    plt.figure(3)
    kok = kav/data[key]['k']
    plt.scatter(R,kok)
    plt.figure(4)
    kvar = ksq - kav*kav
    plt.scatter(R,np.sqrt(kvar)/kav)
    plt.figure(5)
    kdev = np.sqrt(kvar)
    kskew = (kcb - 3*kav*kvar - kav*kav*kav)/(kdev*kdev*kdev)
    plt.scatter(R,kskew)
plt.show()
