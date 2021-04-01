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
        elif row['Unnamed: 1'] == 't':
            thisPHIq = phiq_t(q,row,n)
        elif row['Unnamed: 1']=='aphid':
            thisPHIq = phiq_aphid(q,t,row,n)
        elif row['Unnamed: 1']=='aphi':
            thisPHIq = phiq_aphi(q,row,n)
        elif row['Unnamed: 1'] == 'ax':
            thisPHIq = phiq_ax(q,row,n)
        elif row['Unnamed: 1']=='ay':
            thisPHIq = phiq_ay(q,row,n)
        elif row['Unnamed: 1']=='rt':
            thisPHIq = phiq_rt(q,row,n)
        else:
            print('This constraint type is not defined')
        PHIq=np.vstack((PHIq,thisPHIq))

    return PHIq[1:]
def phiq_ax(q,row,n):
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
    R=np.mat([[0,-1],[1,0]])
    
    thisPHIq = x.T*R*Ai*si
    thisPHIq = np.hstack((x.T,thisPHIq))
    
    if i>0:
        thisPHIql = np.mat([0]*i*3)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*(n-i*3-3))
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_ay(q,row,n):
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
    R=np.mat([[0,-1],[1,0]])
    
    thisPHIq = y.T*R*Ai*si
    thisPHIq = np.hstack((y.T,thisPHIq))
    
    if i>0:
        thisPHIql = np.mat([0]*i*3)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*(n-i*3-3))
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_rt(q,row,n):
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
    Bi=R*Ai
    Aij=Ai.T*Aj

    PHIqi=np.hstack(((-1)*v.T*Bi.T,(-1)*v.T*Ai.T*(rj-ri)-v.T*Aij*sj))/np.linalg.norm(v)
    PHIqj=np.hstack((v.T*Bi.T,v.T*Aij*sj))/np.linalg.norm(v)

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
            thisPHIql = np.mat([0]*j*3)
            thisPHIq = np.hstack((thisPHIql,thisPHIq))

        if i*3+3<n :
            thisPHIqr = np.mat([0]*2*(n-i*3-3))
            thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq


def phiq_t(q,row,n):
    i=row['B_i']
    xi=q[i*3,0]
    yi=q[i*3+1,0]
    phii=q[i*3+2,0]
    ri=np.mat([xi,yi]).T
    Ai=np.mat([[cos(phii),-sin(phii)],[sin(phii),cos(phii)]])
    numbers = row['Unnamed: 6'].split(' ')
    numbers = [ float(x) for x in numbers ]
    si=np.mat(numbers).T
    numbers = row['Unnamed: 7'].split(' ')
    numbers = [ float(x) for x in numbers ]
    vi=np.mat(numbers).T
    
    j=row['B_j']
    xj=q[j*3,0]
    yj=q[j*3+1,0]
    phij=q[j*3+2,0]
    rj=np.mat([xj,yj]).T
    Aj=np.mat([[cos(phij),-sin(phij)],[sin(phij),cos(phij)]])
    numbers = row['Unnamed: 3'].split(' ')
    numbers = [ float(x) for x in numbers ]
    sj=np.mat(numbers).T
    numbers = row['Unnamed: 4'].split(' ')
    numbers = [ float(x) for x in numbers ]
    vj=np.mat(numbers).T

    c=float(row['Unnamed: 8'])
    R=np.mat([[0,-1],[1,0]])
    Bi=R*Ai
    Aij=Ai.T*Aj

    # 算PHIqi
    item1=-1*vi.T*Bi.T
    item2=np.mat([0,0])
    PHIqil=np.vstack((item1,item2))

    item1=-1*vi.T*Ai.T*(rj-ri)-vi.T*Aij*sj
    item2=-1*vi.T*Aij*vj
    PHIqir=np.vstack((item1,item2))

    PHIqi=np.hstack((PHIqil,PHIqir))

    # 算PHIqj
    item1=vi.T*Bi.T
    item2=np.mat([0,0])
    PHIqjl=np.vstack((item1,item2))

    item1=vi.T*Aij*sj
    item2=vi.T*Aij*vj
    PHIqjr=np.vstack((item1,item2))

    PHIqj=np.hstack((PHIqjl,PHIqjr))

    if i<j:
        if i+1==j:
            thisPHIq = np.hstack((PHIqi,PHIqj))
        else:
            thisPHIq = np.hstack((PHIqi,np.mat([0]*(j-i-1)*6).reshape(2,-1),PHIqj))
        
        if i>0:
            thisPHIql = np.mat([0]*i*6).reshape(2,-1)
            thisPHIq = np.hstack((thisPHIql,thisPHIq))

        if j*3+3<n :
            thisPHIqr = np.mat([0]*2*(n-j*3-3)).reshape(2,-1)
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
            thisPHIqr = np.mat([0]*2*(n-j*3-3)).reshape(2,-1)
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
        thisPHIql = np.mat([0]*i*3)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*(n-i*3-3))
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

def phiq_aphi(q,row,n):
    i=row['B_i']
    thisPHIq = np.mat([0,0,1])
    if i>0:
        thisPHIql = np.mat([0]*i*3)
        thisPHIq = np.hstack((thisPHIql,thisPHIq))

    if i*3+3<n :
        thisPHIqr = np.mat([0]*(n-i*3-3))
        thisPHIq = np.hstack((thisPHIq,thisPHIqr))
    return thisPHIq

if __name__ == "__main__":
    q=np.mat([1,0,0,0,0,0]).T

    t=1
    constraints=pd.read_excel('constraints.xlsx')[1:]
    print(calculate_PHIq(q, t, constraints))