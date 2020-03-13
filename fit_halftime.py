import numpy as np
import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def sigmoid(t, A, k, t_half):
    return A/(1 + np.exp(-k*(t - t_half)))

with open('data/vanDerSchootCrowds15_runs100.json') as file:
    data15 = json.load(file)
    moments15 = data15['data']['moments']

t15 = [i['t'] for i in moments15 ]
mass15 = [i['M'] for i in moments15]

init_Akt15 = [mass15[-1], .001, 2500]
best_vals15, covar15 = curve_fit(sigmoid, t15, mass15, p0=init_Akt15)
print("A, k, T = {}".format(best_vals15))

with open('data/vanDerSchootCrowds075.json') as file:
    data075 = json.load(file)
    moments075 = data075['data']['moments']

t075 = [i['t'] for i in moments075 ]
mass075 = [i['M'] for i in moments075]

init_Akt075 = [mass075[-1], .001, 2500]
best_vals075, covar075 = curve_fit(sigmoid, t075, mass075, p0=init_Akt075)
print("A, k, T = {}".format(best_vals075))

with open('data/vanDerSchootCrowds0375.json') as file:
    data0375 = json.load(file)
    moments0375 = data0375['data']['moments']

t0375 = [i['t'] for i in moments0375 ]
mass0375 = [i['M'] for i in moments0375]

init_Akt0375 = [mass0375[-1], .001, 2500]
best_vals0375, covar0375 = curve_fit(sigmoid, t0375, mass0375, p0=init_Akt0375)
print("A, k, T = {}".format(best_vals0375))

with open('data/vanDerSchootNoCrowds.json') as file:
    data0 = json.load(file)
    moments0 = data0['data']['moments']

t0 = [i['t'] for i in moments0 ]
mass0 = [i['M'] for i in moments0]

init_Akt0 = [mass0[-1], .001, 2500]
best_vals0, covar0 = curve_fit(sigmoid, t0, mass0, p0=init_Akt0)
print("A, k, T = {}".format(best_vals0))

sigm0 = [sigmoid(i, best_vals0[0], best_vals0[1], best_vals0[2]) for i in t0]
sigm0375 = [sigmoid(i, best_vals0375[0], best_vals0375[1], best_vals0375[2]) for i in t0375]
sigm075 = [sigmoid(i, best_vals075[0], best_vals075[1], best_vals075[2]) for i in t075]
sigm15 = [sigmoid(i, best_vals15[0], best_vals15[1], best_vals15[2]) for i in t15]

t_half = [np.log(best_vals0[2]), np.log(best_vals0375[2]), np.log(best_vals075[2]), np.log(best_vals15[2])]
phi = [0,0.0375,0.075,0.15]


plt.figure()
plt.scatter(t15, mass15, c='r', marker='+')
plt.scatter(t075, mass075, c='b', marker='+')
plt.scatter(t0375, mass0375, c='g', marker='+')
plt.scatter(t0, mass0, c='k', marker='+')
plt.plot(t15, sigm15, 'r')
plt.plot(t075, mass075, 'b')
plt.plot(t0375, mass0375, 'g')
plt.plot(t0, mass0, 'k')
plt.figure()
plt.plot(phi,t_half)
plt.show()
