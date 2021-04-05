# First developed by JED on March 28, 2021
# As the major assignment for a course

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from calculate_PHI import calculate_PHI
from calculate_PHIq import calculate_PHIq
from calculate_V import calculate_V
from calculate_Gamma import calculate_Gamma
from draw_pic import draw_pic
from show_animation import show_animation
from numpy.linalg import inv,det   #逆矩阵，行列式
from math import *


def main(argv=None):

    #  Read the constraints (make sure all constraints are input normatively)  #
    
    constraints=pd.read_excel('constraints.xlsx')[1:]

    #  *******************Please Input Simulation Parameters*****************  #

    t0 = 0
    N = 300
    te = 6
    maxi = 10000
    e1 = 1e-6
    e2 = 1e-6
    dt = (te-t0)/N

    #返回i的总数乘以3，即约束的个数记录在n中   
    di=constraints['B_i']
    dj=constraints['B_j']
    df=pd.concat([di,dj], ignore_index=True)
    n = df.value_counts().shape[0]*3 

    #  **********************Estimate Initial Condition**********************  #

    # q0 = np.mat(np.random.random(n))
    # 如果上面那个初始化不行，请用户在下面修改
    q0=[0,0,0,0,-2,pi/3,1.5,0,pi/3,0,2,-pi/6]

    #  ******************************Main Program****************************  #

    # 记录每个t下的位形，速度，加速度
    series_q = np.mat([0]*n).T
    series_dq = np.mat([0]*n).T
    series_ddq = np.mat([0]*n).T

    q = np.mat(q0).T

    time_series=np.linspace(t0, te, N+1)
    for t in time_series:
        PHI = calculate_PHI(q, t, constraints)
        DEL = np.sqrt(PHI.T * PHI)

        i = 1

        # 求解t时刻的位形
        while DEL > e1:
            PHIq = calculate_PHIq(q, t, constraints)

            if abs(det(PHIq)) < e2:
                print('Improper initial value:')
                print('abs(det(PHIq)) < e2')
                print(PHIq)
                print(det(PHIq))
                return
            dq=-inv(PHIq)*PHI
            q = q + dq

            PHI = calculate_PHI(q, t, constraints)
            DEL = np.sqrt(PHI.T * PHI)

            if i>=maxi:
                print('Improper initial value2')
                print('i>=maxi')
                return
            i += 1
        
        series_q = np.hstack((series_q,q))

        # 求解t时刻的速度,加速度

        PHIq = calculate_PHIq(q, t, constraints)
        V = calculate_V(q, t, constraints)

        if abs(det(PHIq)) < e2:
            print('Singularity')
            print(PHIq)
            return
        
        dq = inv(PHIq)*V
        series_dq = np.hstack((series_dq,dq))

        Gamma=calculate_Gamma(q,dq,t, constraints)
        ddq=inv(PHIq)*Gamma
        series_ddq = np.hstack((series_ddq,ddq))

        q = q + dq * dt + ddq * dt * dt / 2

    # context是打包好的参数，为了画图，别动
    context = [time_series,series_q,series_dq,series_ddq]


    #  **************************Draw Pictures************************  #

    # n表示第n个物体（从0开始，与输入表格对应）的位形，默认第0个
    # sign是一个字符串,默认为'x'
    # 如要描述位形，有'x'，'y'，'phi'三种选项，
    # 如要描述速度，有'dx'，'dy'，'dphi'三种选项，
    # 如要描述加速度，有'ddx'，'ddy'，'ddphi'三种选项，
    # draw_pic(context,sign='x',n=0):

    draw_pic(context,'x',3) 
    draw_pic(context,'dx',3) 
    draw_pic(context,'ddx',3) 


    #  **************************Show Animation************************  #
    # 如果有画动画的必要，在geometrical_parameters.xlsx中输入所有刚体的参数
    # 1.动画需要输入刚体的几何参数，每一行代表一个刚体，两个端点的坐标目前只支持平面矩形，连体基的x轴为矩形沿着长边的中轴线
    # 2.如有必要手动调整动画显示的坐标轴范围
    show_animation(series_q,dt)

    


if __name__ == "__main__":
    print("Let's Go")
    main()
    print("0 Warnings, 0 Errors")