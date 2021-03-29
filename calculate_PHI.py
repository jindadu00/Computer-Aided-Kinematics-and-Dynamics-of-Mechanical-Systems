import numpy as np
import pandas as pd
from math import *


def calculate_PHI(q, t, constraints):
    PHI = np.mat([0])
    for index, row in constraints.iterrows():
        if row['Unnamed: 1'] == 'ax':
            thisPHI = phi_ax(q,row)
        if row['Unnamed: 1']=='ay':
            thisPHI = phi_ay(q,row)
        if row['Unnamed: 1'] == 'r':
            thisPHI = phi_r(q,row)
        if row['Unnamed: 1']=='aphid':
            thisPHI = phi_aphid(q,t,row)
        if row['Unnamed: 1'] == 'rt':
            thisPHI = phi_rt(q,row)

        PHI=np.vstack((PHI,thisPHI))
    return PHI[1:]

def phi_ax(q,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    ri=np.mat([xi,yi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    ci=float(row['Unnamed: 8'])
    x=np.mat([1,0]).T
    return x.T*(ri+Ai*si)-ci

def phi_ay(q,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    ri=np.mat([xi,yi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    ci=float(row['Unnamed: 8'])
    y=np.mat([0,1]).T
    return y.T*(ri+Ai*si)-ci


def phi_r(q,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    ri=np.mat([xi,yi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T

    j=row['B_j']
    xj=q[j*3,0]
    yj=q[j*3+1,0]
    phij=q[j*3+2,0]
    rj=np.mat([xj,yj]).T
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T

    return rj+Aj*sj-ri-Ai*si

def phi_rt(q,row):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    ri=np.mat([xi,yi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    
    j=row['B_j']
    xj=q[j*3,0]
    yj=q[j*3+1,0]
    phij=q[j*3+2,0]
    rj=np.mat([xj,yj]).T
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T

    numbers = row['Unnamed: 7'].split(' ')
    numbers = [ float(x) for x in numbers ]
    v=np.mat(numbers).T
    
    c=float(row['Unnamed: 8'])
    R=np.mat([[0,-1],[1,0]])

    h=rj+Aj*sj-ri-Ai*si
    
    return (R*Ai*v).T*h/np.linalg.norm(v)-c

def phi_aphid(q,t,row):
    i=row['B_i']
    phii=q[i*3+2,0]
    c=str(row['Unnamed: 8']).replace("pi", "np.pi").replace("^", "**")
    return np.mat(phii-eval(c))


if __name__ == "__main__":
    q=np.mat([1,0,0,0,0,0]).T

    t=1
    constraints=pd.read_excel('constraints.xlsx')[1:]
    print(calculate_PHI(q, t, constraints))