import numpy as np
import pandas as pd
from math import *

def calculate_V(q, t, constraints):
    V = np.mat([0])
    for index, row in constraints.iterrows():
        if row['Unnamed: 1'] == 'r':
            thisV = v_r(q,row)
        elif row['Unnamed: 1'] == 't':
            thisV = v_t(q,row)
        elif row['Unnamed: 1'] == 'ax':
            thisV = v_ax(q,row)
        elif row['Unnamed: 1'] == 'ay':
            thisV = v_ay(q,row)
        elif row['Unnamed: 1'] == 'rt':
            thisV = v_rt(q,row)
        elif row['Unnamed: 1']=='aphi':
            thisV = v_aphi(q,row)
        elif row['Unnamed: 1']=='aphid':
            thisV = v_aphid(q,t,row)
        else:
            print('This constraint type is not defined')
        V=np.vstack((V,thisV))

    return V[1:]

def v_ax(q,row):
    return np.mat([0]).T
def v_ay(q,row):
    return np.mat([0]).T
def v_aphi(q,row):
    return np.mat([0]).T
def v_rt(q,row):
    return np.mat([0]).T
def v_r(q,row):
    return np.mat([0,0]).T
def v_t(q,row):
    return np.mat([0,0]).T
def v_aphid(q,t,row):
    i=row['B_i']
    c=str(row['Unnamed: 9']).replace("pi", "np.pi").replace("^", "**")
    return np.mat(eval(c))
