import numpy as np
import json
import matplotlib.pyplot as plt

masses = {}
numbers = {}
run_masses = {}
run_numbers = {}
run_ts = {}

with open('data/fluc_halftime_volume/Smol_N10k_1000runs.json') as file:
    data = json.load(file)
    masses['10k'] = {'mass':[i['M']/10000.0 for i in data['data']['moments']]}
    masses['10k']['dev'] = [i['M_dev']/10000.0 for i in data['data']['moments']]
    run_masses['10k'] = [[i['M']/10000.0 for i in run] for run in data['data']['runMoments']]
    numbers['10k'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['10k']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['10k'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_ts['10k'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
    t = [i['t'] for i in data['data']['moments']]
with open('data/fluc_halftime_volume/Smol_N7k_1000runs.json') as file:
    data = json.load(file)
    masses['7k'] = {'mass':[i['M']/7000.0 for i in data['data']['moments']]}
    masses['7k']['dev'] = [i['M_dev']/7000.0 for i in data['data']['moments']]
    numbers['7k'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['7k']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['7k'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['7k'] = [[i['M']/7000.0 for i in run] for run in data['data']['runMoments']]
    run_ts['7k'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N5k_1000runs.json') as file:
    data = json.load(file)
    masses['5k'] = {'mass':[i['M']/5000.0 for i in data['data']['moments']]}
    masses['5k']['dev'] = [i['M_dev']/5000.0 for i in data['data']['moments']]
    numbers['5k'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['5k']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['5k'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['5k'] = [[i['M']/5000.0 for i in run] for run in data['data']['runMoments']]
    run_ts['5k'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N2.5k_1000runs.json') as file:
    data = json.load(file)
    masses['2.5k'] = {'mass':[i['M']/2500.0 for i in data['data']['moments']]}
    masses['2.5k']['dev'] = [i['M_dev']/2500.0 for i in data['data']['moments']]
    numbers['2.5k'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['2.5k']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['2.5k'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['2.5k'] = [[i['M']/2500.0 for i in run] for run in data['data']['runMoments']]
    run_ts['2.5k'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N1k_1000runs.json') as file:
    data = json.load(file)
    masses['1k'] = {'mass':[i['M']/1000.0 for i in data['data']['moments']]}
    masses['1k']['dev'] = [i['M_dev']/1000.0 for i in data['data']['moments']]
    numbers['1k'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['1k']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['1k'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['1k'] = [[i['M']/1000.0 for i in run] for run in data['data']['runMoments']]
    run_ts['1k'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N750_1000runs.json') as file:
    data = json.load(file)
    masses['750'] = {'mass':[i['M']/750.0 for i in data['data']['moments']]}
    masses['750']['dev'] = [i['M_dev']/750.0 for i in data['data']['moments']]
    numbers['750'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['750']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['750'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['750'] = [[i['M']/750.0 for i in run] for run in data['data']['runMoments']]
    run_ts['750'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N500_1000runs.json') as file:
    data = json.load(file)
    masses['500'] = {'mass':[i['M']/500.0 for i in data['data']['moments']]}
    masses['500']['dev'] = [i['M_dev']/500.0 for i in data['data']['moments']]
    numbers['500'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['500']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['500'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['500'] = [[i['M']/500.0 for i in run] for run in data['data']['runMoments']]
    run_ts['500'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N250_1000runs.json') as file:
    data = json.load(file)
    masses['250'] = {'mass':[i['M']/250.0 for i in data['data']['moments']]}
    masses['250']['dev'] = [i['M_dev']/250.0 for i in data['data']['moments']]
    numbers['250'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['250']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['250'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['250'] = [[i['M']/250.0 for i in run] for run in data['data']['runMoments']]
    run_ts['250'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N100_1000runs.json') as file:
    data = json.load(file)
    masses['100'] = {'mass':[i['M']/100.0 for i in data['data']['moments']]}
    masses['100']['dev'] = [i['M_dev']/100.0 for i in data['data']['moments']]
    numbers['100'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['100']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['100'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['100'] = [[i['M']/100.0 for i in run] for run in data['data']['runMoments']]
    run_ts['100'] = [[i['t'] for i in run] for run in data['data']['runMoments']]
with open('data/fluc_halftime_volume/Smol_N50_1000runs.json') as file:
    data = json.load(file)
    masses['50'] = {'mass':[i['M']/50.0 for i in data['data']['moments']]}
    masses['50']['dev'] = [i['M_dev']/50.0 for i in data['data']['moments']]
    numbers['50'] = {'number':[i['P'] for i in data['data']['moments']]}
    numbers['50']['dev'] = [i['P_dev'] for i in data['data']['moments']]
    run_numbers['50'] = [[i['P'] for i in run] for run in data['data']['runMoments']]
    run_masses['50'] = [[i['M']/50.0 for i in run] for run in data['data']['runMoments']]
    run_ts['50'] = [[i['t'] for i in run] for run in data['data']['runMoments']]

plt.figure()
plt.plot(t, masses['100']['mass'], 'r')
for idx, run in enumerate(run_masses['100']):
    plt.plot(run_ts['100'][idx], run)
plt.figure()
for key in masses:
    plt.plot(t, masses[key]['dev'], label=key)
plt.legend()

fig, axs = plt.subplots(3,3)
axs[0, 0].plot(t, masses['10k']['mass'])
axs[0,0].set_title('10k')
for idx, run in enumerate(run_masses['10k']):
    axs[0, 0].plot(run_ts['10k'][idx], run)

axs[0, 1].plot(t, masses['7k']['mass'])
axs[0,1].set_title('7k')
for idx, run in enumerate(run_masses['7k']):
    axs[0, 1].plot(run_ts['7k'][idx], run)

axs[0,2].plot(t, masses['5k']['mass'])
axs[0,2].set_title('5k')
for idx, run in enumerate(run_masses['5k']):
    axs[0,2].plot(run_ts['5k'][idx], run)

axs[1,0].plot(t, masses['2.5k']['mass'])
axs[1,0].set_title('2.5k')
for idx, run in enumerate(run_masses['2.5k']):
    axs[1,0].plot(run_ts['2.5k'][idx], run)

axs[1,1].plot(t, masses['1k']['mass'])
axs[1,1].set_title('1k')
for idx, run in enumerate(run_masses['1k']):
    axs[1,1].plot(run_ts['1k'][idx], run)

axs[1,2].plot(t, masses['750']['mass'])
axs[1,2].set_title('750')
for idx, run in enumerate(run_masses['750']):
    axs[1,2].plot(run_ts['750'][idx], run)

axs[2,0].plot(t, masses['500']['mass'])
axs[2,0].set_title('500')
for idx, run in enumerate(run_masses['500']):
    axs[2,0].plot(run_ts['500'][idx], run)

axs[2,1].plot(t, masses['250']['mass'])
axs[2,1].set_title('250')
for idx, run in enumerate(run_masses['250']):
    axs[2,1].plot(run_ts['250'][idx], run)

axs[2,2].plot(t, masses['100']['mass'])
axs[2,2].set_title('100')
for idx, run in enumerate(run_masses['100']):
    axs[2,2].plot(run_ts['100'][idx], run)

fig2, axs2 = plt.subplots(3,3)
axs2[0, 0].plot(t, numbers['10k']['number'])
axs2[0,0].set_title('10k')
for idx, run in enumerate(run_numbers['10k']):
    axs2[0, 0].plot(run_ts['10k'][idx], run)

axs2[0, 1].plot(t, numbers['7k']['number'])
axs2[0,1].set_title('7k')
for idx, run in enumerate(run_numbers['7k']):
    axs2[0, 1].plot(run_ts['7k'][idx], run)

axs2[0,2].plot(t, numbers['5k']['number'])
axs2[0,2].set_title('5k')
for idx, run in enumerate(run_numbers['5k']):
    axs2[0,2].plot(run_ts['5k'][idx], run)

axs2[1,0].plot(t, numbers['2.5k']['number'])
axs2[1,0].set_title('2.5k')
for idx, run in enumerate(run_numbers['2.5k']):
    axs2[1,0].plot(run_ts['2.5k'][idx], run)

axs2[1,1].plot(t, numbers['1k']['number'])
axs2[1,1].set_title('1k')
for idx, run in enumerate(run_numbers['1k']):
    axs2[1,1].plot(run_ts['1k'][idx], run)

axs2[1,2].plot(t, numbers['750']['number'])
axs2[1,2].set_title('750')
for idx, run in enumerate(run_numbers['750']):
    axs2[1,2].plot(run_ts['750'][idx], run)

axs2[2,0].plot(t, numbers['500']['number'])
axs2[2,0].set_title('500')
for idx, run in enumerate(run_numbers['500']):
    axs2[2,0].plot(run_ts['500'][idx], run)

axs2[2,1].plot(t, numbers['250']['number'])
axs2[2,1].set_title('250')
for idx, run in enumerate(run_numbers['250']):
    axs2[2,1].plot(run_ts['250'][idx], run)

axs2[2,2].plot(t, numbers['100']['number'])
axs2[2,2].set_title('100')
for idx, run in enumerate(run_numbers['100']):
    axs2[2,2].plot(run_ts['100'][idx], run)

plt.show()