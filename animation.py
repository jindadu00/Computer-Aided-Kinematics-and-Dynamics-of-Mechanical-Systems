# coding:UTF-8

'''
Kinematics analysis
Build by Luo Xueling

输入文件示例：INP2.txt
_______________________________________________________________________
Number of rigidbodys = 5
0	2	0	1.5
0	0	0	4
2	2	0	0.2
0	4	0	1.9
0	4	0	0.2

Number of constraints = 11
ax		0	0 0		0 0							1 0
ay		0	0 0		0 0							1 2
ax		1	0 0		0 0							1 0
ay		1	0 0		0 0							1 0
r		0	1.5	0	0 0		2	0 0		0 0
t		1	0 0		1 0		2	0 0		1 0
r		1	4 0		0 0		3	1.9 0	0 0
r		3	0 0		0 0		4	0 0		0 0
ay		4	0 0		0 0							1 4
aphi	4	0 0		0 0							1 0
aphid	0	0 0		0 0							2 0.44 3.1415926

t0 = 0
te = 6
dt = 0.001
e1 = 0.0001
e2 = 0.0001
_______________________________________________________________________
输入文件说明：
对于刚体，其输入为
-----------------------
x0	y0	phi0	length
-----------------------
其中，length用于作图，初始坐标全为0可能导致雅可比矩阵奇异，应当输入估计值
对于约束，其输入为
---------------------------------------------------------------------------------------------------------------------------------------------------
type(lower case)	body1	s1'T	v1'T	body2		s2'T	v2'T	c(t)项数	c(t)各项系数（幂次由小到大）	C（ad,add约束可能用到的坐标点）
---------------------------------------------------------------------------------------------------------------------------------------------------
Number of constraints不指约束方程个数
'''

import os
import math
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import random
import shutil
import imageio,imageio_ffmpeg
import numpy as np
import multiprocessing
from functools import partial
import sys

plt.rcParams['figure.dpi'] = 300


def run():
    # bgProcess程序由c++编写，编译为
    bgProcess()
    outputPNG(inpFilePath)
    outputMP4(inpFilePath)
    return "Python ends"


def bgProcess():
    # https://blog.csdn.net/weixin_33946605/article/details/85139653
    # https://cloud.tencent.com/developer/ask/115790
    selfpath = os.path.realpath(sys.argv[0]).replace("plot.py", "")
    exe_path = selfpath + "Kinematics.exe"
    params = "\"" + root + " " + inpname + "\\\""
    with open('./runBgProcess.bat', 'w') as file:
        # 不生成新
        file.write("@echo off\n")
        file.write("powershell.exe -command ^\n")
        file.write("  \"start " + exe_path + " \\" + params + "\" -NoNewWindow -Wait")
    os.system('runBgProcess.bat')


# 做曲线图
def outputPNG(inpFilePath):
    print("Python function (outputPNG) start")
    print("Reading " + inpFilePath)
    # 坐标轴刻度设置 https://www.jb51.net/article/164187.htm
    # df.plot https://blog.csdn.net/h_hxx/article/details/90635650
    # 子图的标题、标签设置 https://zhuanlan.zhihu.com/p/55675217
    # 输入文件路径预处理
    bodyNum = getBodyNum()
    fig = plt.figure(figsize=(15, 4 * bodyNum))
    for i in range(bodyNum):
        pospath = root + inpname + "_OUT_body" + str(i) + "_pos.csv"
        velpath = root + inpname + "_OUT_body" + str(i) + "_vel.csv"
        accpath = root + inpname + "_OUT_body" + str(i) + "_acc.csv"
        path = {'Position': pospath, 'Velocity': velpath, 'Accelaration': accpath}
        labelname = []
        labelname.append(['x', 'y', 'phi', 'x,y/(m)', 'phi/(rad)'])
        labelname.append(['dx/dt', 'dy/dt', 'dphi/dt', 'Vx,Vy/(m/s)', 'Vphi/(rad/s)'])
        labelname.append(['ddx/ddt', 'ddy/ddt', 'ddphi/ddt', 'Ax,Ay/(m/s2)', 'Aphi/(rad/s2)'])
        keys = list(path.keys())

        for cnt in range(3):
            df = pd.read_csv(path[keys[cnt]], header=None, names=['t', 'x', 'y', 'phi'])
            df = df.set_index(df['t'].values)

            ax1 = fig.add_subplot(bodyNum, 3, 3 * i + cnt + 1)
            df['x'].plot(ax=ax1, label=labelname[cnt][0], color='Blue')
            df['y'].plot(ax=ax1, label=labelname[cnt][1], color='Red')
            ax1.set_ylabel(labelname[cnt][3])
            plt.legend(loc="upper left")

            ax2 = ax1.twinx()
            df['phi'].plot(ax=ax2, label=labelname[cnt][2], color='Green')
            ax2.set_ylabel(labelname[cnt][4])

            ymax = max(df['x'].max(), df['y'].max())
            ymin = min(df['x'].min(), df['y'].min())
            m = max(abs(ymax), abs(ymin))
            if m != 0:
                ax1.set_ylim([-1.5 * m, 1.5 * m])

            m = max(abs(df['phi'].max()), abs(df['phi'].min()))
            if m != 0:
                ax2.set_ylim([-1.5 * m, 1.5 * m])

            plt.legend(loc="upper right")
            plt.title('Body ' + str(i) + "- " + keys[cnt] + " plot")
            plt.xlabel('t/s')

    plt.tight_layout()
    fig.savefig(root + inpname + "_OUT_Plot.png", dpi=300)
    print('outputPNG completed')


# MP4主程序
def outputMP4(inpFilePath):
    # https://blog.csdn.net/nkhgl/article/details/103185573
    # https://blog.csdn.net/qq_28888837/article/details/85778395
    # https://vimsky.com/zh-tw/examples/detail/python-method-imageio.mimsave.html
    # https://www.cnpython.com/qa/458636

    print("Python function (outputMP4) start")
    print("Reading " + inpFilePath)
    bodyNum = getBodyNum()
    folder = root + inpname + "_Anim//"

    # 每个刚体长度，及相关点坐标
    length = []
    with open(root + inpname + '.txt') as file:
        file.readline()
        for i in range(bodyNum):
            line = file.readline()
            st = line.split('\t')[-1].split(' ')[-1]
            st = st.replace("\n", "")
            length.append(float(st))

    # 读取数据文件
    df = []

    for i in range(bodyNum):
        pospath = root + inpname + "_OUT_body" + str(i) + "_pos.csv"
        dfi = pd.read_csv(pospath, header=None, names=['t', 'x', 'y', 'phi'])
        dfi = dfi.set_index(dfi['t'].values)

        df.append(dfi)

    # 获取各时间点坐标信息
    xmin = -1
    ymin = -1
    xmax = 1
    ymax = 1
    points = []
    totaltime = df[0]['t'].values[-1]
    n = len(df[0]['t'].values)
    for t in range(n):
        tmp = [[], []]

        for body in range(bodyNum):
            x1 = df[body]['x'].values[t]
            y1 = df[body]['y'].values[t]
            phi = df[body]['phi'].values[t]
            x2 = x1 + length[body] * math.cos(phi)
            y2 = y1 + length[body] * math.sin(phi)
            tmp[0].append(x1)
            tmp[0].append(x2)
            tmp[1].append(y1)
            tmp[1].append(y2)
        if min(tmp[0]) < xmin: xmin = min(tmp[0])
        if min(tmp[1]) < ymin: ymin = min(tmp[1])
        if max(tmp[0]) > xmax: xmax = max(tmp[0])
        if max(tmp[1]) > ymax: ymax = max(tmp[1])
        points.append(tmp)

    # 随机生成刚体颜色
    def randomClr():
        randclr = (random.random(), random.random(), random.random())
        return randclr

    clr = []
    for body in range(bodyNum): clr.append(randomClr())

    # 保持axis equal
    d = 1.3333 * (ymax - ymin) / (xmax - xmin)
    if d > 1:
        xl = (xmin) * 1.3 * d
        xr = (xmax) * 1.3 * d
        yl = (ymin) * 1.3
        yr = (ymax) * 1.3
    else:
        xl = (xmin) * 1.3
        xr = (xmax) * 1.3
        yl = (ymin) * 1.3 / d
        yr = (ymax) * 1.3 / d

    if os.path.exists(folder):
        shutil.rmtree(folder)

    print("Plotting image.")
    os.makedirs(folder)

    # https://www.itranslater.com/qa/details/2582631824283403264
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=math.floor(cores))
    tlist = np.linspace(0, n - 1, n, dtype=int)
    poolableFunc = partial(threadPng, df, bodyNum, folder, points, n, xl, xr, yl, yr, clr)
    pool.map(poolableFunc, tlist)

    # 输出mp4
    print("Saving mp4")
    frames = []
    for t in range(n):
        frames.append(imageio.imread(folder + str(t) + ".png"))
    # imageio.mimsave(root + inpname + '_output.gif', frames, 'GIF', duration=totaltime/10/n)
    imageio.mimsave(root + inpname + '_output.mp4', frames, 'MP4', fps=n / totaltime)
    # shutil.rmtree(folder)

    print('outputMP4 completed')


# MP4中生成Png的函数，用于多线程运行
def threadPng(df, bodyNum, folder, points, n, xl, xr, yl, yr, clr, t):
    fig, ax = plt.subplots()
    currenttime = df[0]['t'].values[t]
    ax.set_xlim([xl, xr])
    ax.set_ylim([yl, yr])
    ax.grid('on')
    # 作刚体
    for body in range(bodyNum):
        xy0 = points[t][0][2 * body:2 * body + 2]
        xy1 = points[t][1][2 * body:2 * body + 2]

        line = mlines.Line2D(xy0, xy1, lw=(xr - xl) * 0.7, ls='-', color=clr[body])
        line.set_zorder(0)
        ax.add_line(line)

    # 作铰
    x = points[t][0]
    y = points[t][1]
    for p in range(len(x)):
        circle = mpatches.Circle((x[p], y[p]), radius=(xr - xl) * 0.005 * math.log(bodyNum + 5), ec="blue", fc='white')
        circle.set_zorder(1)
        ax.add_patch(circle)

    ax.set_title("Time = " + str(round(currenttime, 5)) + " s")
    ax.set_xlim([xl, xr])
    ax.set_ylim([yl, yr])

    # 保存并记录图像
    fig.savefig(folder + str(t) + ".png", dpi=300)
    print("Image " + str(len(os.listdir(folder))) + "/" + str(n) + " has been plotted")
    fig.clf()
    plt.close(fig)


# 路径预处理
def pathPreProcess(inpFilePath):
    # 字符串预处理
    st = inpFilePath.split("/")
    root = ""
    for part in range(len(st) - 1):
        root += (st[part])
        root += (r"/")
    # 含后缀名
    filename = st[-1]
    # 不含后缀名
    inpname = filename.split(".")[0]
    return root, filename, inpname


def getBodyNum():
    dirlist = list(os.walk(root))[0][2]
    # 输出文件名
    outnamelist = []

    for item in dirlist:
        if (inpname in item.split('_')) and (".txt" not in item) and (".png" not in item):
            outnamelist.append(item)

    bodyNum = math.floor(len(outnamelist) / 3)

    return bodyNum


if __name__ == "__main__":
    print("Python start")
    print("Input the path of INP file.")
    print("---------------Example---------------")
    print("G://KinematicsExer//INP1.txt")
    print("-------------------------------------")
    inpFilePath = input()
    root, filename, inpname = pathPreProcess(inpFilePath)
    bodyNum = getBodyNum()
    run()
