import glob, os, json
import pandas as pd
import numpy as np

def reactionData(jsondata):
    colnames = [label for label in jsondata['data']['reactions'][0] if label !='dt' if label != 't']
    r_data = pd.DataFrame(jsondata['data']['reactions'], columns=colnames)
    t_data = pd.DataFrame(jsondata['data']['reactions'], columns=['t'])
    r_data = r_data.div(r_data.sum(axis=1), axis=0)
    return r_data.join(t_data['t'])

def runData(jsondata):
    outData={}
    for idx,run in enumerate(jsondata['data']['runMoments']):
        outData[idx]={}
        outData[idx]['M'] = [i['M'] for i in run]
        outData[idx]['P'] = [i['P'] for i in run]
        outData[idx]['L'] = [i['L'] for i in run]
        outData[idx]['t'] = [i['t'] for i in run]
    return outData

def getData(rel_folder_path, data_key, runs=False):
    data={}
    files = glob.glob(rel_folder_path+'/*.json')
    for filepath in files:
        with open(filepath) as file:
            jsondata=json.load(file)
            try:
                key = str(jsondata[data_key])
            except:
                key = '0'
            N = jsondata['N']
            data[key] = {}
            try:
                data[key]['a'] = jsondata['a']
            except:
                pass
            try:
                data[key]['b'] = jsondata['b']
            except:
                pass
            try:
                data[key]['ka'] = jsondata['ka']
            except:
                pass
            try:
                data[key]['kb'] = jsondata['kb']
            except:
                pass
            try:
                data[key]['kn'] = jsondata['kn']
            except:
                pass
            try:
                data[key]['k2'] = jsondata['k2']
            except:
                pass
            try:
                data[key]['nc'] = jsondata['nc']
            except:
                pass
            try:
                data[key]['n2'] = jsondata['n2']
            except:
                pass
            data[key]['N'] = jsondata['N']
            data[key]['ind_runs']=jsondata['ind_runs']
            data[key]['conc']=jsondata['Co']
            data[key]['M']=[i['M'] for i in jsondata['data']['moments']]
            data[key]['M_norm']=[i/N for i in data[key]['M']]
            data[key]['M2']=[i['M2'] for i in jsondata['data']['moments']]
            data[key]['M_dev']=[i['M_dev']/N for i in jsondata['data']['moments']]
            data[key]['M_cv']=[i['M_dev']/i['M'] if i['M']!=0 else 0 for i in jsondata['data']['moments']]
            data[key]['P']=[i['P'] for i in jsondata['data']['moments']]
            data[key]['P2']=[i['P2'] for i in jsondata['data']['moments']]
            data[key]['P_dev']=[i['P_dev'] for i in jsondata['data']['moments']]
            data[key]['P_cv']=[i['P_dev']/i['P'] if i['P']!=0 else 0 for i in jsondata['data']['moments']]
            data[key]['L']=[i['L'] for i in jsondata['data']['moments']]
            data[key]['L2']=[i['L2'] for i in jsondata['data']['moments']]
            data[key]['L_dev']=[i['L_dev'] for i in jsondata['data']['moments']]
            data[key]['L_cv']=[i['L_dev']/i['L'] if i['L']!=0 else 0 for i in jsondata['data']['moments']]
            data[key]['t']=[i['t'] for i in jsondata['data']['moments']]
            colnames = [label for label in jsondata['data']['reactions'][0] if label !='dt' if label != 't']
            r_data = pd.DataFrame(jsondata['data']['reactions'], columns=colnames)
            t_data = pd.DataFrame(jsondata['data']['reactions'], columns=['t'])
            sum_rows = pd.DataFrame()
            sum_rows['sum'] = r_data.sum(axis=1) 
            r_data = r_data.div(r_data.sum(axis=1), axis=0)
            r_data = r_data.join(t_data['t'])
            r_data = r_data.join(sum_rows['sum'])
            data[key]['reactions'] = r_data
            data[key]['histograms'] = jsondata['data']['histograms']
            if runs:
                data[key]['runs'] = runData(jsondata)

    return data

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
def sigmoid(t, A, k, t_half):
    output = A/(1 + np.exp(-k*(t - t_half)))
    return output

def half_index(mass):
    half_val = mass[-1]/2.0
    for idx, m in enumerate(mass):
        if m >= half_val:
            return idx
    return 0.1