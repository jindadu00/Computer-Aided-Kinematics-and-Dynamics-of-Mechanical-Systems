import numpy as np
import pandas as pd
from math import *

def calculate_PHIq(q, t, constraints):
    PHIq = np.mat([0]*q.shape[0])

    di=constraints['B_i']
    dj=constraints['B_j']
    df=pd.concat([di,dj], ignore_index=True)
    n = df.value_counts().shape[0]*3

    for index, row in constraints.iterrows():
        if row['Unnamed: 1'] == 'r':
            thisPHIq = phiq_r(q,row,n)

        if row['Unnamed: 1']=='aphid':
            thisPHIq = phiq_aphid(q,t,row,n)

        if row['Unnamed: 1'] == 'ax':
            thisPHIq = phiq_ax(q,row,n)

        if row['Unnamed: 1']=='ay':
            thisPHIq = phiq_ay(q,row,n)
            
        PHIq=np.vstack((PHIq,thisPHIq))

    return PHIq[1:]

def phiq_ax(q,row,n):
    i=row['B_i']
    phii=q[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    
    x=np.mat([1,0]).T
    R=np.mat([[0,-1],[1,0]])
    
    thisPHIq = x.T*R*Ai*si
    thisPHIq = np.hstack((x.T,thisPHIq))
    
    if i>0:
        thisPHIql = np.mat([0]*i*6).reshape(2,-1)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*2*(n-i*3-3)).reshape(1,-1)
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_ay(q,row,n):
    i=row['B_i']
    phii=q[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    
    y=np.mat([0,1]).T
    R=np.mat([[0,-1],[1,0]])
    
    thisPHIq = y.T*R*Ai*si
    thisPHIq = np.hstack((y.T,thisPHIq))
    
    if i>0:
        thisPHIql = np.mat([0]*i*6).reshape(2,-1)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*2*(n-i*3-3)).reshape(1,-1)
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_r(q,row,n):
    i=row['B_i']
    phii=q[i*3+2,0]
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T

    j=row['B_j']
    phij=q[j*3+2,0]
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T
    
    R=np.mat([[0,-1],[1,0]])
    PHIqi=-1*(np.hstack((np.identity(2),R*Ai*si)))
    PHIqj=np.hstack((np.identity(2),R*Aj*sj))

    if i<j:
        if i+1==j:
            thisPHIq = np.hstack((PHIqi,PHIqj))
        else:
            thisPHIq = np.hstack((PHIqi,np.mat([0]*(j-i-1)*6).reshape(2,-1),PHIqj))
        
        if i>0:
            thisPHIql = np.mat([0]*i*6).reshape(2,-1)
            thisPHIq = np.hstack((thisPHIql,thisPHIq))

        if j*3+3<n :
            thisPHIqr = np.mat([0]*2*(n-j*3-3)).reshape(1,-1)
            thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    else:
        if j+1==i:
            thisPHIq = np.hstack((PHIqj,PHIqi))
        else:
            thisPHIq = np.hstack((PHIqj,np.mat([0]*(i-j-1)*6).reshape(2,-1),PHIqi))
        
        if j>0:
            thisPHIql = np.mat([0]*j*6).reshape(2,-1)
            thisPHIq = np.hstack((thisPHIql,thisPHIq))

        if i*3+3<n :
            thisPHIqr = np.mat([0]*2*(n-i*3-3)).reshape(1,-1)
            thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_aphid(q,t,row,n):
    i=row['B_i']
    thisPHIq = np.mat([0,0,1])
    if i>0:
        thisPHIql = np.mat([0]*i*6).reshape(2,-1)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*2*(n-i*3-3)).reshape(1,-1)
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

if __name__ == "__main__":
    q=np.mat([1,0,0,0,0,0]).T

    t=1
    constraints=pd.read_excel('constraints.xlsx')[1:]
    print(calculate_PHIq(q, t, constraints))