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
models = ['secondary_nuc','2d', 'frag']

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
    
    # plt.plot(crange, cube_spl(crange))
    poly_cs = np.polyfit(cs, ts, 1)
    poly = np.poly1d(poly_cs)
    poly_slope = poly.deriv()
    plt.plot(crange, poly(crange))
    # l_1 = np.polyfit(c_sorted[:9],t_sorted[:9],1)
    # line_1 = np.poly1d(l_1)
    # line_1_slope = line_1.deriv()
    # s1 = line_1_slope(0)
    # plt.plot([-.2,1.63], line_1([-.2,1.63]), 'r--')
    # l_2 = np.polyfit(c_sorted[16:20],t_sorted[16:20],1)
    # line_2 = np.poly1d(l_2)
    # line_2_slope = line_2.deriv()
    # s2 = line_2_slope(0)
    # plt.plot([1.63,3.75], line_2([1.63,3.75]), 'r--')
    plt.scatter(cs,ts)
    # l_3 = np.polyfit(c_sorted[-4:],t_sorted[-4:],1)
    # line_3 = np.poly1d(l_3)
    # line_3_slope = line_3.deriv()
    # s3 = line_3_slope(0)
    # plt.plot([3.75, 4.9], line_3([3.75,4.9]), 'r--')
    # plt.text(.35,2,'$\gamma$ = {}'.format(round(s1,2)))
    # plt.text(2.8, -1.7, '$\gamma$ = {}'.format(round(s2,2)))
    # plt.text(4.2,-4.4, '$\gamma$ = {}'.format(round(s3,2)))
    plt.scatter(cs,ts, color='blue')
    #plt.title('$ln(t_{1/2})$ vs $ln(c_0)$')
    # plt.axvline(x=1.2, color='k', linestyle = '--')
    # plt.axvline(x=2.4, color='k', linestyle = '--')
    # plt.text(.25, -5.8, 'Region 1', weight='bold')
    # plt.text(1.4,-5.8, 'Region 2', weight='bold')
    # plt.text(2.55, -5.8, 'Region 3', weight='bold')
    plt.ylabel('$ln$ $t_{1/2}$')
    plt.xlabel('$ln$ $c_0$')
    plt.title(folder)
    plt.figure()
    plt.plot(c_sorted, poly_slope(c_sorted))
    plt.title(folder)
    # plt.figure(2)
    # plt.plot(crange, poly_slope(crange))

    # plt.figure(1)
    # t_slope_unsort = [(t-ts[idx])/(cs[idx+1]-cs[idx]) for idx,t in enumerate(ts[1:])]
    # c_avg = [(c+cs[idx])/2 for idx, c in enumerate(cs[1:])]
    # t_slope = [t for _,t in sorted(zip(c_avg, t_slope_unsort))]
    # c_avg.sort()
    # logc = np.arange(c_avg[0], c_avg[-1], (c_avg[-1]-c_avg[0])/100)
    # plt.plot(c_avg, t_slope)
    # plt.plot()

# data['reactions'] = {}

# data['reactions']['0.9um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_0.9.12add025-5b1c-4f71-86ae-0bf01041e8bf.json')
# data['reactions']['1um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_1.b3f4f491-a370-4203-a085-7d2d92045af0.json')
# data['reactions']['2um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_2_4 (1).1c879311-5ceb-47b1-95ea-702a7f0fece9.json')

# fig, axes = plt.subplots(nrows=3, ncols=3)



# data['reactions']['2um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[0,2], legend=False)
# data['reactions']['1um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[0,1], legend=False)
# data['reactions']['0.9um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[0,0], legend=False)
# axes[0,2].set_xlim(0,4)
# #ax4.set_title('relative reaction rate vs t ($c_0$=4um)')
# axes[0,1].set_xlim(0,8)
# axes[0,1].set_title('Relative reaction rate vs t')
# axes[0,0].set_xlim(0,20)
# axes[0,0].set_xticks([0,10,20])
# axes[0,0].set_yticks([0,.2,.4,.6,.8])
# axes[0,1].set_xticks([0,4,8])
# axes[0,1].set_yticks([0,.2,.4,.6,.8])
# axes[0,2].set_xticks([0,2,4])
# axes[0,2].set_yticks([0,.2,.4,.6,.8])
# axes[0,0].text(10,0.7,'$c_0$ = 0.9$\mu$M', size='10')
# axes[0,1].text(4,0.7,'$c_0$ = 1$\mu$M', size='10')
# axes[0,2].text(2,0.7,'$c_0$ = 2$\mu$M', size='10')
# #ax1.set_title('relative reaction rate vs t ($c_0$=1um)')
# #ax4.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))



# data['reactions']['3um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_3.0604ef80-a769-4f8c-a860-f665312ecd14.json')
# data['reactions']['4um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_2_4 (2).3d7868f6-8daf-497b-b755-0cea50ceb61d.json')
# data['reactions']['6.5um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_6.5.059f3b21-2365-401f-850a-2dcca7324e8b.json')


# data['reactions']['3um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[1,0], legend=False)
# data['reactions']['4um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[1,1], legend=False)
# data['reactions']['6.5um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[1,2], legend=False)
# #ax14.set_title('relative reaction rate vs t ($c_0$=14um)')
# #ax18.set_title('relative reaction rate vs t ($c_0$=18um)')
# #ax30.set_title('relative reaction rate vs t ($c_0$=30um)')
# axes[1,0].set_xticks([0,6,12])
# axes[1,0].set_yticks([0,.2,.4,.6,.8])
# axes[1,1].set_xticks([0,6,12])
# axes[1,1].set_yticks([0,.2,.4,.6,.8])
# #axes[1,2].set_xticks([0,.2,.4])
# axes[1,2].set_yticks([0,.2,.4,.6,.8])
# axes[1,0].text(6,0.7,'$c_0$ = 3$\mu$M', size='10')
# axes[1,1].text(6,0.7,'$c_0$ = 4$\mu$M', size='10')
# axes[1,2].text(1.2,0.75,'$c_0$ = 6.5$\mu$M', size='10')


# data['reactions']['30um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_30.2b21e245-af57-45c2-b3b8-8781e9709639.json')
# data['reactions']['60um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_60.bb5d99d2-7508-47aa-baa8-8727aa2bb513.json')
# data['reactions']['100um'] = reactionData('data/scaling_laws/secondary_nuc/2dary_nuc_c0_100.ad965326-216d-4f30-bd60-c45a426e7e60.json')
# data['reactions']['30um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[2,0], legend=False)
# data['reactions']['60um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[2,1], legend=False)
# data['reactions']['100um'].plot(x='t',y=['nucleation','addition','fragmentation','secondary-nucleation'], ax=axes[2,2], legend=False)

# # ax80.set_title('relative reaction rate vs t ($c_0$=80um)')
# # ax100.set_title('relative reaction rate vs t ($c_0$=100um)')
# # ax120.set_title('relative reaction rate vs t ($c_0$=120um)')
# #axes[2,0].set_xticks([0,.03,.06,])
# axes[2,0].set_yticks([0,.2,.4,.6,.8])
# #axes[2,1].set_xticks([0,.02,.04, .06])
# axes[2,1].set_yticks([0,.2,.4,.6,.8])
# #axes[2,2].set_xticks([0,.01,.02,.03])
# axes[2,2].set_yticks([0,.2,.4,.6,.8])
# axes[2,0].text(.19,0.78,'$c_0$ = 30$\mu$M', size='10')
# axes[2,1].text(.04,0.84,'$c_0$ = 60$\mu$M', size='10')
# axes[2,2].text(.025,0.86,'$c_0$ = 100$\mu$M', size='10')

# for ax in axes[0,:]:
#     axis = ax.axes.get_xaxis()
#     label = axis.get_label()
#     label.set_visible(False)
# for ax in axes[1,:]:
#     axis = ax.axes.get_xaxis()
#     label = axis.get_label()
#     label.set_visible(False)
# for ax in axes[2,:]:
#     axis = ax.axes.get_xaxis()
#     label = axis.get_label()
#     label.set_visible(False)
# axes[2,1].axes.get_xaxis().get_label().set_visible(True)
# axes[2,1].set_label('t (s)')
# fig.tight_layout()
# fig.subplots_adjust(hspace=0.2, wspace=0.25)

plt.show()

