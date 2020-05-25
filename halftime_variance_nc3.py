import numpy as np
import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import glob

def sigmoid(t, A, k, t_half):
    output = A/(1 + np.exp(-k*(t - t_half)))
    return output

def half_index(mass):
    half_val = mass[-1]/2.0
    for idx, m in enumerate(mass):
        if m >= half_val:
            return idx
    return 0.1
def half_time_line(ms,ts, half_mass):
    m = (ms[1]-ms[0])/(ts[1]-ts[0])
    b = ms[0]-m*ts[0]
    t_half = (half_mass - b)/m
    return t_half

def simple_half_time(M, t, half_mass):
    for idx,m in enumerate(M):
        if m > half_mass:
            return t[idx]

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

# masses = {}
# with open('data/fluc_halftime_volume/Smol_N10k_1000runs.json') as file:
#     data = json.load(file)
#     masses['10k']=[i['M'] for i in data['data']['moments']]
#     t = [i['t'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N7k_1000runs.json') as file:
#     data = json.load(file)
#     masses['7k']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N5k_1000runs.json') as file:
#     data = json.load(file)
#     masses['5k']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N2.5k_1000runs.json') as file:
#     data = json.load(file)
#     masses['2.5k']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N1k_1000runs.json') as file:
#     data = json.load(file)
#     masses['1k']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N750_1000runs.json') as file:
#     data = json.load(file)
#     masses['750']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N500_1000runs.json') as file:
#     data = json.load(file)
#     masses['500']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N250_1000runs.json') as file:
#     data = json.load(file)
#     masses['250']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N100_1000runs.json') as file:
#     data = json.load(file)
#     masses['100']=[i['M'] for i in data['data']['moments']]
# with open('data/fluc_halftime_volume/Smol_N50_1000runs.json') as file:
#     data = json.load(file)
#     masses['50']=[i['M'] for i in data['data']['moments']]

ddir = 'data/fluc_halftime_volume/nc3'
files = glob.glob(ddir+'/*.json')
data={}
for filepath in files:
    with open(filepath) as file:
        print(filepath)
        jsondata = json.load(file)
        key = str(jsondata['N'])
        data[key] = {}
        data[key]['ind_runs']=jsondata['ind_runs']
        data[key]['conc']=jsondata['Co']
        data[key]['mass']=[i['M'] for i in jsondata['data']['moments']]
        t=[i['t'] for i in jsondata['data']['moments']]
        data[key]['t']=t
        data[key]['t_half']=half_time(data[key]['mass'],t)
        data[key]['N_inv']=1/jsondata['N']
        data[key]['N']=jsondata['N']
        mMax = data[key]['mass'][-1]
        half_guess = t[half_index(data[key]['mass'])]
        k_guess = mMax/(40*half_guess)
        guess = [mMax, k_guess, half_guess]
        data[key]['fit'],_ = curve_fit(sigmoid, t, data[key]['mass'], guess)
        data[key]['runs']={}
        data[key]['th_sq']=0
        data[key]['th_av']=0
        for idx,run in enumerate(jsondata['data']['runMoments']):
            data[key]['runs'][idx]=[]
            data[key]['runs'][idx].append([step['t'] for step in run])
            data[key]['runs'][idx].append([step['M'] for step in run])
            # M_half = mMax/2
            # htime = simple_half_time(data[key]['runs'][idx][1],data[key]['runs'][idx][0], M_half)
            # mMax = data[key]['runs'][idx][1][-1]
            half_guess = data[key]['runs'][idx][0][half_index(data[key]['runs'][idx][1])]
            k_guess = 23*mMax/(half_guess)
            guess = [mMax, k_guess, half_guess]
            fit,_ = curve_fit(sigmoid, t, data[key]['mass'], guess)
            data[key]['th_sq']+=fit[2]**2
            data[key]['th_av']+=fit[2]
        data[key]['th_sq'] = data[key]['th_sq']/data[key]['ind_runs']
        data[key]['th_av'] = (data[key]['th_av']/data[key]['ind_runs'])**2
        data[key]['th_var'] = data[key]['th_sq'] - data[key]['th_av']
        data[key]['th_dev'] = np.sqrt(data[key]['th_var'])
        print(data[key]['th_var'])

th_var_multiple = [data[key]['th_var'] for key in data if data[key]['N'] % 3 == 0]
th_var_1 = [data[key]['th_var'] for key in data if data[key]['N'] % 3 == 1]
th_var_2 = [data[key]['th_var'] for key in data if data[key]['N'] % 3 == 2]
th_dev_multiple = [data[key]['th_dev'] for key in data if data[key]['N'] % 3 == 0]
th_dev_1 = [data[key]['th_dev'] for key in data if data[key]['N'] % 3 == 1]
th_dev_2 = [data[key]['th_dev'] for key in data if data[key]['N'] % 3 == 2]
th_CV_multiple = [data[key]['th_dev']/data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 0]
th_CV_1 = [data[key]['th_dev']/data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 1]
th_CV_2 = [data[key]['th_dev']/data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 2]
N_inv_multiple = [data[key]['N_inv'] for key in data if data[key]['N'] % 3 == 0]
N_inv_1 = [data[key]['N_inv'] for key in data if data[key]['N'] % 3 == 1]
N_inv_2 = [data[key]['N_inv'] for key in data if data[key]['N'] % 3 == 2]
t_half_multiple = [data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 0]
t_half_1 = [data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 1]
t_half_2 = [data[key]['fit'][2] for key in data if data[key]['N'] % 3 == 2]
t_multiple_sorted = [t for _,t in sorted(zip(N_inv_multiple,t_half_multiple))]
t_1_sorted = [t for _,t in sorted(zip(N_inv_1,t_half_1))]
t_2_sorted = [t for _,t in sorted(zip(N_inv_2,t_half_2))]
tvar_multiple_sorted = [t for _,t in sorted(zip(N_inv_multiple,th_var_multiple))]
tvar_1_sorted = [t for _,t in sorted(zip(N_inv_1,th_var_1))]
tvar_2_sorted = [t for _,t in sorted(zip(N_inv_2,th_var_2))]
N_inv_multiple_sorted = sorted(N_inv_multiple)
N_inv_1_sorted = sorted(N_inv_1)
N_inv_2_sorted = sorted(N_inv_2)
l_1e = np.polyfit(N_inv_multiple_sorted,t_multiple_sorted,1)
line_1e = np.poly1d(l_1e)
l_1o = np.polyfit(N_inv_1_sorted,t_1_sorted,1)
line_1o = np.poly1d(l_1o)
l_1x = np.polyfit(N_inv_2_sorted,t_2_sorted,1)
line_1x = np.poly1d(l_1x)
l_2e = np.polyfit(N_inv_multiple_sorted, tvar_multiple_sorted,1)
line_2e = np.poly1d(l_2e)
l_2o = np.polyfit(N_inv_1_sorted, tvar_1_sorted,1)
line_2o = np.poly1d(l_2o)
l_2x = np.polyfit(N_inv_2_sorted, tvar_2_sorted,1)
line_2x = np.poly1d(l_2x)
plt.figure()
plt.scatter(N_inv_multiple_sorted, t_multiple_sorted, label='0')
plt.scatter(N_inv_1_sorted, t_1_sorted, label='1')
plt.scatter(N_inv_2_sorted, t_2_sorted, label='2')
#plt.xlim(left=0, right=.021)
# plt.ylim(bottom=-.11)
plt.title('$t_{1/2}$ vs 1/N')
plt.legend()
plt.plot(N_inv_multiple_sorted, line_1e(N_inv_multiple_sorted), 'r--', linewidth=0.5)
plt.plot(N_inv_1_sorted, line_1o(N_inv_1_sorted), 'r--', linewidth=0.5)
plt.plot(N_inv_2_sorted, line_1x(N_inv_2_sorted), 'r--', linewidth=0.5)
plt.figure()
plt.scatter(N_inv_multiple, th_var_multiple, label='0')
plt.scatter(N_inv_1,th_var_1, label='1')
plt.scatter(N_inv_2,th_var_2, label='2')
plt.legend()
# plt.plot(N_inv_sorted, line_2(N_inv_sorted), 'r--', linewidth=0.5)
#plt.xlim(left=0,right=0.0205)
plt.ylim(bottom=0)
plt.figure()
plt.scatter(N_inv_multiple, th_dev_multiple, label='0')
plt.scatter(N_inv_1,th_dev_1, label='1')
plt.scatter(N_inv_2,th_dev_2, label='2')
plt.legend()
#plt.xlim(left=0,right=0.0205)
plt.ylim(bottom=0)
plt.figure()
plt.scatter(N_inv_multiple,th_CV_multiple, label='0')
plt.scatter(N_inv_1,th_CV_1, label='1')
plt.scatter(N_inv_2,th_CV_2, label='2')
plt.legend()
#plt.xlim(left=0,right=0.0205)
plt.ylim(bottom=0)

# for key, dat in data.items():
#     plt.figure()
#     plt.plot(dat['t'],dat['mass'], linewidth=2)
#     for idx, run in dat['runs'].items():
#         plt.plot(run[0],run[1],linestyle='-.', linewidth=0.5)
#         plt.title('N = {}'.format(key))
# plt.show()
# fits = {}
# fit_params = {}
# N = [10000,7000,5000,2500,1000,750,500,250,100,50]
# N_inv = [1/i for i in N]


# t_halves = [data[key]['t_half'] for key in data]
# log_t_halves = [np.log(i) for i in t_halves]
# plt.figure()
# plt.scatter(N_inv,log_t_halves)
plt.show()