from matplotlib import pyplot as plt

#解决中文显示问题
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def draw_pic(context,sign='x',n=0):
    # n表示第n个物体（从0开始，与输入表格对应）的位形，默认第0个
    # sign是一个字符串,默认为'x'
    # 如要描述位形，有'x'，'y'，'phi'三种选项，
    # 如要描述速度，有'dx'，'dy'，'dphi'三种选项，
    # 如要描述加速度，有'ddx'，'ddy'，'ddphi'三种选项，
    if sign=='x' or sign=='y' or sign=='phi':
        draw_q(context[1],context[0],sign,n)
    if sign=='dx' or sign=='dy' or sign=='dphi':
        draw_dq(context[2],context[0],sign,n)
    if sign=='ddx' or sign=='ddy' or sign=='ddphi':
        draw_ddq(context[3],context[0],sign,n)

def draw_q(series_q,time_series,sign,n):

    if sign=='x':
        i = 0
    if sign=='y':
        i = 1
    if sign=='phi':
        i = 2

    t = time_series.reshape(-1,1)
    q = series_q[i+n*3,1:].T
    
    plt.grid()
    plt.title("第"+str(n)+"个刚体的"+sign+"-t曲线")
    plt.xlabel('t(s)')
    if sign=='phi':
        sign=sign+'(rad)'
    else:
        sign=sign+'(cm)'
    plt.ylabel(sign)
    plt.plot(t,q)
    plt.show()

def draw_dq(series_dq,time_series,sign,n):

    if sign=='dx':
        i = 0
    if sign=='dy':
        i = 1
    if sign=='dphi':
        i = 2
    t = time_series.reshape(-1,1)
    dq = series_dq[i+n*3,1:].T
    plt.grid()
    plt.title("第"+str(n)+"个刚体的"+sign+"-t曲线")
    plt.xlabel('t(s)') 
    if sign=='dphi':
        sign=sign+'(rad/s)'
    else:
        sign=sign+'(cm/s)'
    plt.ylabel(sign)

    plt.plot(t,dq)
    plt.show()

def draw_ddq(series_ddq,time_series,sign,n):

    if sign=='ddx':
        i = 0
    if sign=='ddy':
        i = 1
    if sign=='ddphi':
        i = 2
    t = time_series.reshape(-1,1)
    ddq = series_ddq[i+n*3,1:].T
    plt.grid()
    plt.title("第"+str(n)+"个刚体的"+sign+"-t曲线")
    plt.xlabel('t(s)') 
    if sign=='ddphi':
        sign=sign+'(rad/s^2)'
    else:
        sign=sign+'(cm/s^2)'
    plt.ylabel(sign)
    plt.plot(t,ddq)
    plt.show()