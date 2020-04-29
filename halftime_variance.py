import numpy as np
import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def sigmoid(t, A, k, t_half):
    return A/(1 + np.exp(-k*(t - t_half)))

def half_index(mass):
    half_val = mass[-1]/2.0
    for idx, m in enumerate(mass):
        print(m)
        if m >= half_val:
            return idx
    return 0.1

masses = {}
with open('data/fluc_halftime_volume/Smol_N10k_1000runs.json') as file:
    data = json.load(file)
    masses['10k']=[i['M'] for i in data['data']['moments']]
    t = [i['t'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N7k_1000runs.json') as file:
    data = json.load(file)
    masses['7k']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N5k_1000runs.json') as file:
    data = json.load(file)
    masses['5k']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N2.5k_1000runs.json') as file:
    data = json.load(file)
    masses['2.5k']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N1k_1000runs.json') as file:
    data = json.load(file)
    masses['1k']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N750_1000runs.json') as file:
    data = json.load(file)
    masses['750']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N500_1000runs.json') as file:
    data = json.load(file)
    masses['500']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N250_1000runs.json') as file:
    data = json.load(file)
    masses['250']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N100_1000runs.json') as file:
    data = json.load(file)
    masses['100']=[i['M'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N50_1000runs.json') as file:
    data = json.load(file)
    masses['50']=[i['M'] for i in data['data']['moments']]

fits = {}
fit_params = {}
N = [10000,7000,5000,2500,1000,750,500,250,100,50]
N_inv = [1/i for i in N]

for key in masses:
    A = masses[key][-1]
    t_half_guess = half_index(masses[key])
    k_guess = 2.0/t_half_guess
    init_Akt = [A, k_guess, t_half_guess]
    params, covar = curve_fit(sigmoid, t, masses[key], p0=init_Akt)
    fit_params[key] = params
    fits[key] = [sigmoid(i, params[0], params[1], params[2]) for i in t]


t_halves = [fit_params[key][2] for key in fit_params]
log_t_halves = [np.log(i) for i in t_halves]
plt.figure()
plt.scatter(N_inv,log_t_halves)
plt.show()