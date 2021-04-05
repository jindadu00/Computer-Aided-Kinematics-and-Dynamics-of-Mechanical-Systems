from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import matplotlib.lines as mlines
import pandas as pd


def show_animation(series_q,dt):
    q = series_q[:,1:]

    geometrical_parameters=pd.read_excel('geometrical_parameters.xlsx')

    fig = plt.figure()

    #如有必要手动调整坐标轴范围
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-10,10), ylim=(-10, 10))
    ax.set_aspect('equal')
    ax.grid()

    color=['b','g','r','c','m','y','k']
    lines = [mlines.Line2D([], [], ls='-',color=color[x]) for x in range(geometrical_parameters.shape[0])]

    for line in lines:
        ax.add_line(line)

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def animate(i):
        for k, row in geometrical_parameters.iterrows():
            line=lines[k]
            thisx = [q[k*3,i]+row['xl']*cos(q[2+k*3,i]), q[k*3,i]+row['xr']*cos(q[2+k*3,i])]
            thisy = [q[1+k*3,i]+row['xl']*sin(q[2+k*3,i]), q[1+k*3,i]+row['xr']*sin(q[2+k*3,i])]
            line.set_data(thisx, thisy)
        return lines
    # 通过调整dt*3000的系数可以改变快慢，改变range(100)的参数可以改变帧数
    ani = animation.FuncAnimation(fig, animate, range(100),
                                interval=dt*3000, blit=True, init_func=init)
    ani.save("test.gif",writer='pillow')
    plt.show()