import numpy as np
import json
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import glob, os
import pandas as pd
import matplotlib.gridspec as gridspec


def sigmoid(t, A, k, t_half):
    return A/(1 + np.exp(-k*(t - t_half)))

def half_index(mass):
    half_val = mass[-1]/2.0
    for idx, m in enumerate(mass):
        print(m)
        if m >= half_val:
            return idx
    return 0.1

def half_time_line(ms,ts, half_mass):
    m = (ms[1]-ms[0])/(ts[1]-ts[0])
    b = ms[0]-m*ts[0]
    t_half = (half_mass - b)/m
    return t_half

def half_time(M, t):
    M_half = M[-1]/2.0
    M_4 = 2*M[-1]/5.0
    M_6 = 3*M[-1]/5.0
    idx = 0
    m = M[idx]
    while m < M_4:
        idx = idx+1
        m = M[idx]
    idx_4 = idx
    while m < M_6:
        idx = idx+1
        m = M[idx]
    idx_6 = idx
    return half_time_line([M_4,M_6],[t[idx_4],t[idx_6]],M_half)



def reactionData(path):
    with open(path) as file:
        jsondata = json.load(file)
        colnames = [label for label in jsondata['data']['reactions'][0] if label !='dt' if label != 't']
        r_data = pd.DataFrame(jsondata['data']['reactions'], columns=colnames)
        t_data = pd.DataFrame(jsondata['data']['reactions'], columns=['t'])
        r_data = r_data.div(r_data.sum(axis=1), axis=0)
    return r_data.join(t_data['t'])

data_directory=os.getcwd()+'/data/scaling_laws/'
models = ['secondary_nuc','2d', 'frag','cr_0','cr_0.1','cr_0.2','cr_0.15','comp','comp2']

data = {}

for folder in models:
    data[folder] = {}
    ddir = data_directory+folder
    files = glob.glob(ddir+'/*.json')
    for filepath in files:
        with open(filepath) as file:
            jsondata = json.load(file)
            key = str(jsondata['Co'])
            data[folder][key] = {}
            data[folder][key]['conc']=jsondata['Co']
            data[folder][key]['mass']=[i['M'] for i in jsondata['data']['moments']]
            
            t=[i['t'] for i in jsondata['data']['moments']]
            data[folder][key]['t']=t
            data[folder][key]['t_half']=half_time(data[folder][key]['mass'],t)



# plt.figure(2)
t_slope = {}
for folder in models:
    if folder=='cr_0':
        fig,ax = plt.subplots(1)
        cs = []
        ts = []
        for key in data[folder]:
            cs.append(np.log(data[folder][key]['conc']))
            ts.append(np.log(data[folder][key]['t_half']))
        t_sorted = [t for _,t in sorted(zip(cs,ts))]
        c_sorted = sorted(cs)
        # cube_spl = CubicSpline(c_sorted, t_sorted)
        crange=np.arange(c_sorted[0], c_sorted[-1], (c_sorted[-1]-c_sorted[0])/100)
        poly_cs = np.polyfit(cs, ts, 1)
        poly = np.poly1d(poly_cs)
        poly_slope = poly.deriv()
        ax.plot(crange, poly(crange))
        ax.scatter(cs,ts, color='blue')
        ax.set_ylabel('$ln$ $t_{1/2}$')
        ax.set_xlabel('$ln$ $c_0$')
        phi = folder.split('_')[1]
        ax.text(4, -.8 , '$\phi$ = {}'.format(phi))
        # ax.set_title(folder)
        # axes[0].figure()
        # axes[0].plot(c_sorted, poly_slope(c_sorted))
        # axes[0].title(folder)
    elif folder=='cr_0.1':
        cs = []
        ts = []
        for key in data[folder]:
            cs.append(np.log(data[folder][key]['conc']))
            ts.append(np.log(data[folder][key]['t_half']))
        t_sorted = [t for _,t in sorted(zip(cs,ts))]
        c_sorted = sorted(cs)
        # cube_spl = CubicSpline(c_sorted, t_sorted)
        crange=np.arange(c_sorted[0], c_sorted[-1], (c_sorted[-1]-c_sorted[0])/100)
        poly_cs = np.polyfit(cs, ts, 2)
        poly = np.poly1d(poly_cs)
        poly_slope = poly.deriv()
        ax.plot(crange, poly(crange))
        ax.scatter(cs,ts, color='red')
        phi = folder.split('_')[1]
        ax.text(4, -2 , '$\phi$ = {}'.format(phi))
        # ax.set_ylabel('$ln$ $t_{1/2}$')
        # ax.set_xlabel('$ln$ $c_0$')
        # ax.set_title(folder)
        # axes[1].figure()
        # axes[1].plot(c_sorted, poly_slope(c_sorted))
        # axes[1].title(folder)
    elif folder == 'cr_0.15':
        cs = []
        ts = []
        for key in data[folder]:
            cs.append(np.log(data[folder][key]['conc']))
            ts.append(np.log(data[folder][key]['t_half']))
        t_sorted = [t for _,t in sorted(zip(cs,ts))]
        c_sorted = sorted(cs)
        # cube_spl = CubicSpline(c_sorted, t_sorted)
        crange=np.arange(c_sorted[0], c_sorted[-1], (c_sorted[-1]-c_sorted[0])/100)
        poly_cs = np.polyfit(cs, ts, 2)
        poly = np.poly1d(poly_cs)
        poly_slope = poly.deriv()
        ax.plot(crange, poly(crange),'k')
        ax.scatter(cs,ts, color='black')
        phi = folder.split('_')[1]
        ax.text(4, -3.2 , '$\phi$ = {}'.format(phi))
        # ax.set_ylabel('$ln$ $t_{1/2}$')
        # ax.set_xlabel('$ln$ $c_0$')
        # ax.set_title(folder)
        # axes[1].figure()
        # axes[1].plot(c_sorted, poly_slope(c_sorted))
        # axes[1].title(folder)
    elif folder == 'cr_0.2':
        cs = []
        ts = []
        for key in data[folder]:
            cs.append(np.log(data[folder][key]['conc']))
            ts.append(np.log(data[folder][key]['t_half']))
        t_sorted = [t for _,t in sorted(zip(cs,ts))]
        c_sorted = sorted(cs)
        # cube_spl = CubicSpline(c_sorted, t_sorted)
        crange=np.arange(c_sorted[0], c_sorted[-1], (c_sorted[-1]-c_sorted[0])/100)
        poly_cs = np.polyfit(cs, ts, 1)
        poly = np.poly1d(poly_cs)
        poly_slope = poly.deriv()
        ax.plot(crange, poly(crange))
        ax.scatter(cs,ts, color='green')
        phi = folder.split('_')[1]
        ax.text(4, -4.7 , '$\phi$ = {}'.format(phi))
        # ax.set_ylabel('$ln$ $t_{1/2}$')
        # ax.set_xlabel('$ln$ $c_0$')
        # ax.set_title(folder)
        # axes[1].figure()
        # axes[1].plot(c_sorted, poly_slope(c_sorted))
        # axes[1].title(folder)
    else:       
        plt.figure()
        plt.title(folder)
        cs = []
        ts = []
        for key in data[folder]:
            cs.append(np.log(data[folder][key]['conc']))
            ts.append(np.log(data[folder][key]['t_half']))
        t_sorted = [t for _,t in sorted(zip(cs,ts))]
        c_sorted = sorted(cs)
        cube_spl = CubicSpline(c_sorted, t_sorted)
        crange=np.arange(c_sorted[0], c_sorted[-1], (c_sorted[-1]-c_sorted[0])/100)
        if folder == 'comp2':
            poly_cs = np.polyfit(cs, ts, 2)
        elif folder == 'comp':
            poly_cs = np.polyfit(cs, ts, 2)
        else:
            poly_cs = np.polyfit(cs, ts, 1)
        poly = np.poly1d(poly_cs)
        poly_slope = poly.deriv()
        plt.plot(crange, poly(crange))
        plt.scatter(cs,ts)
        plt.scatter(cs,ts, color='blue')
        plt.ylabel('$ln$ $t_{1/2}$')
        plt.xlabel('$ln$ $c_0$')
        plt.title(folder)
        plt.figure()
        plt.plot(c_sorted, poly_slope(c_sorted))
        plt.title(folder)

plt.show()

