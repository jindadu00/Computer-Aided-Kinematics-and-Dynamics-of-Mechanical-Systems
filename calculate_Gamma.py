import numpy as np
import pandas as pd
from math import *

def calculate_Gamma(q, dq, t, constraints):
    Gamma = np.mat([0])
    for index, row in constraints.iterrows():
        if row['Unnamed: 1'] == 'r':
            thisGamma = gamma_r(q,dq,row)
        elif row['Unnamed: 1'] == 't':
            thisGamma = gamma_t(q,dq,row)
        elif row['Unnamed: 1']=='aphid':
            thisGamma = gamma_aphid(q,t,row)

        elif row['Unnamed: 1']=='ax':
            thisGamma = gamma_ax(q,dq,row)

        elif row['Unnamed: 1']=='ay':
            thisGamma = gamma_ay(q,dq,row)

        elif row['Unnamed: 1']=='rt':
            thisGamma = gamma_rt(q,dq,row)
        else:
            print('This constraint type is not defined')
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

def gamma_t(q,dq,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    dxi=dq[i*3,0]
    dyi=dq[i*3+1,0]
    dphii=dq[i*3+2,0]
    ri=np.mat([xi,yi]).T
    dri=np.mat([dxi,dyi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])

    numbers = row['Unnamed: 7'].split(' ')
    numbers = [ float(x) for x in numbers ]
    vi=np.mat(numbers).T
    
    j=row['B_j']
    xj=q[j*3,0]
    yj=q[j*3+1,0]
    phij=q[j*3+2,0]
    dxj=dq[j*3,0]
    dyj=dq[j*3+1,0]

    rj=np.mat([xj,yj]).T
    drj=np.mat([dxj,dyj]).T
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]

    c=float(row['Unnamed: 8'])
    R=np.mat([[0,-1],[1,0]])
    Bi=R*Ai
    item1=dphii*dphii*vi.T*Bi.T*(rj-ri)+2*dphii*vi.T*Ai.T*(drj-dri)
    item2=np.mat([0])
    return np.vstack((item1,item2))



def gamma_rt(q,dq,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    dxi=dq[i*3,0]
    dyi=dq[i*3+1,0]
    dphii=dq[i*3+2,0]
    ri=np.mat([xi,yi]).T
    dri=np.mat([dxi,dyi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    
    j=row['B_j']
    xj=q[j*3,0]
    yj=q[j*3+1,0]
    phij=q[j*3+2,0]
    dxj=dq[j*3,0]
    dyj=dq[j*3+1,0]
    dphij=dq[j*3+2,0]
    rj=np.mat([xj,yj]).T
    drj=np.mat([dxj,dyj]).T
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T

    numbers = row['Unnamed: 7'].split(' ')
    numbers = [ float(x) for x in numbers ]
    v=np.mat(numbers).T
    c=float(row['Unnamed: 8'])
    R=np.mat([[0,-1],[1,0]])
    Bi=R*Ai
    Bij=R*Ai.T*Aj
    item1=(-2)*v.T*Ai.T*(drj-dri)*dphii
    item2=(-1)*v.T*Bi.T*(rj-ri)*dphii*dphii
    item3=v.T*Bij*sj*(dphij*dphij-dphii*dphii)*(dphij*dphij-dphii*dphii)
    return (item1+item2+item3)*(-1/np.linalg.norm(v))

def gamma_aphi(q,row):
    return np.mat([0])

def gamma_aphid(q,t,row):
    c=str(row['Unnamed: 10']).replace("pi", "np.pi").replace("^", "**")
    return np.mat(eval(c))