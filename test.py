# # First developed by JED on March 28, 2021
# # As the major assignment for a course

# import numpy as np
# import pandas as pd
# from matplotlib import pyplot as plt
# from calculate_PHI import calculate_PHI
# from calculate_PHIq import calculate_PHIq
# from calculate_V import calculate_V
# from calculate_Gamma import calculate_Gamma
# from draw_pic import draw_pic
# from numpy.linalg import inv,det   #逆矩阵，行列式
# from math import *


# def main(argv=None):
#     constraints=pd.read_excel('constraints.xlsx')[1:]

#     t0 = 0
#     N = 1
#     te = 0.5
#     maxi = 100000
#     e1 = 1e-4
#     e2 = 1e-4
#     dt = (te-t0)/N

#     #返回i的总数乘以3，即约束的个数记录在n中   
#     di=constraints['B_i']
#     dj=constraints['B_j']
#     df=pd.concat([di,dj], ignore_index=True)
#     n = df.value_counts().shape[0]*3 

#     q0 = np.mat(np.random.random(n))
#     q = np.mat(q0).T
#     dq = q/2
#     time_series=np.linspace(t0, te, N+1)
#     t=t0
#     PHI = calculate_PHI(q, t, constraints)
#     PHIq = calculate_PHIq(q, t, constraints)
#     V = calculate_V(q, t, constraints)
#     Gamma=calculate_Gamma(q,dq,t, constraints)
#     # print('PHI')
#     # print(PHI.shape)
#     # print('PHIq')
#     # print(PHIq.shape)
#     # print('V')
#     # print(V.shape)
#     # print('Gamma')
#     # print(Gamma.shape)

# if __name__ == "__main__":
#     main()
print('hello\\nworld')