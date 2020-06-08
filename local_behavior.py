import numpy as np
import json
import matplotlib.pyplot as plt
from utilities import getData
from matplotlib.ticker import FormatStrFormatter
import pandas as pd

def find_idx(data, val):
    for idx,i in enumerate(data):
        if i>val:
            return idx-1

conc_folder = 'data/crowders-local/conc'
crowd_folder = 'data/crowders-local'

concData = getData(conc_folder, 'Co', runs=True)
crowdData = getData(crowd_folder, 'phi', runs=True)

fig, axes = plt.subplots(1,2)
pidx=0
for key,data in concData.items():
    N=data['N']
    axes[pidx].plot(data['t'][:30],[i/N for i in data['M'][:30]])
    for tkey,run in data['runs'].items():
        idx = find_idx(run['t'], data['t'][30])
        axes[pidx].plot(run['t'][:idx],[i/N for i in run['M'][:idx]], linestyle='--', linewidth=0.5)
    axes[pidx].set_title('M vs t ($c_0$ = {})'.format(key))
    axes[pidx].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axes[pidx].set_ylim(top=1)
    pidx += 1

fig2, axes2 = plt.subplots(2,2)
xidx=0
yidx=0
for key,data in crowdData.items():
    N = data['N']
    axes2[yidx,xidx].plot(data['t'],[i/N for i in data['M']])
    for tkey,run in data['runs'].items():
        # idx = find_idx(run['t'], data['t'])
        if int(tkey)%4==0:
            axes2[yidx,xidx].plot(run['t'],[i/N for i in run['M']], linestyle='--', linewidth=0.5)
    axes2[yidx,xidx].set_title('M vs t ($\phi$ = {})'.format(key))
    axes2[yidx,xidx].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axes2[yidx,xidx].set_ylim(top=1)
    xidx += 1
    if xidx % 2 == 0:
        xidx = 0
        yidx += 1

fig3, axes3 = plt.subplots(2,2)
xidx=0
yidx=0
leg = False
lines = ['nucleation','addition','fragmentation','secondary-nucleation']
markers = ['o','s','x','^']
for key,data in crowdData.items():
    if 2*yidx + xidx == 1:
        leg = True
    else:
        leg = False
    for reac,marker in zip(lines,markers):
        axes3[yidx,xidx].scatter(data['reactions']['t'],data['reactions'][reac], marker=marker, s=2)
    if leg:
        axes3[yidx,xidx].legend(loc=1, prop={'size': 5})
        axes3[yidx,xidx].text(0.85, 0.65, '$\phi$ = {}'.format(key),
            verticalalignment='top', horizontalalignment='center',
            transform=axes3[yidx,xidx].transAxes, fontsize=8)
    else:
        axes3[yidx,xidx].text(0.85, 0.95, '$\phi$ = {}'.format(key),
            verticalalignment='top', horizontalalignment='center',
            transform=axes3[yidx,xidx].transAxes, fontsize=8)

    xidx += 1
    if xidx % 2 == 0:
        xidx=0
        yidx+=1
plt.show()