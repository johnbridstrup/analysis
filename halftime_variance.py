import numpy as np
import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def sigmoid(t, A, k, t_half):
    return A/(1 + np.exp(-k*(t - t_half)))

with open('data/VanDerSchoot_phi0_50ind_N1000.json') as file:
    data = json.load(file)
    moments = data['data']['runMoments']

masses = []
ts = []
fits = []
fit_params = []

for moment in moments:
    run_mass = [i['M'] for i in moment]
    run_t = [i['t'] for i in moment]
    ts.append(run_t)
    masses.append(run_mass)
    half_idx = 0
    half_val = run_mass[-1] / 2.0;
    for idx, mom in enumerate(moment):
        if mom['M'] > half_val:
            half_idx = idx
            break
    
    init_Akt = [run_mass[-1], .001, run_t[half_idx]]
    params, covar = curve_fit(sigmoid, run_t, run_mass, p0=init_Akt)
    fit_params.append(params)
    fits.append([sigmoid(i, params[0], params[1], params[2]) for i in run_t])

t_spread = [i[2] for i in fit_params]
plt.hist(t_spread)
plt.show()