import numpy as np
import pandas as pd
from math import *

def calculate_Gamma(q, dq, t, constraints):
    Gamma = np.mat([0])
    for index, row in constraints.iterrows():
        if row['Unnamed: 1'] == 'r':
            thisGamma = gamma_r(q,dq,row)
            
        if row['Unnamed: 1']=='aphid':
            thisGamma = gamma_aphid(q,t,row)

        if row['Unnamed: 1']=='ax':
            thisGamma = gamma_ax(q,dq,row)

        if row['Unnamed: 1']=='ay':
            thisGamma = gamma_ay(q,dq,row)

        Gamma=np.vstack((Gamma,thisGamma))

    return Gamma[1:]

def gamma_ax(q,dq,row):
    i=row['B_i']
    phii=q[i*3+2,0]
    dphii=dq[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    x=np.mat([1,0]).T

    return x.T*Ai*si*dphii*dphii

def gamma_ay(q,dq,row):
    i=row['B_i']
    phii=q[i*3+2,0]
    dphii=dq[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    y=np.mat([0,1]).T

    return y.T*Ai*si*dphii*dphii

def gamma_r(q,dq,row):
    i=row['B_i']
    phii=q[i*3+2,0]
    dphii=dq[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T

    j=row['B_j']
    phij=q[j*3+2,0]
    dphij=dq[j*3+2,0]
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T

    return Aj*sj*dphij*dphij-Ai*si*dphii*dphii

def gamma_aphid(q,t,row):
    c=str(row['Unnamed: 10']).replace("pi", "np.pi").replace("^", "**")
    return np.mat(eval(c))