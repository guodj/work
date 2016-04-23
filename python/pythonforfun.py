#!/home/gdj/anaconda3/bin/python3
#-*- coding: utf-8 -*-

__author__ = 'Guo Dongjie'

"""
Some funny and useful python programmes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def f1():
    """
    复利的力量，假设每月存x元，年利率为p，问
    n年后翻多少倍？
    """
    x = 2000
    p = 0.2
    n = 10

    s = 0
    for k1 in range(n):
        temp = x*((1+p)**(n-k1))
        for k2 in range(12):
            temp = temp - x*p/12*k2
            s = s + temp
    print('每月存{}元,共存{}年,设年利率为{:.0f}%'.format(x,n,p*100))
    print("则总存款数为{},"
          "共收回{:>.2f}元,是存款的{:>.1f}倍。".format(x*12*n,s,s/(x*n*12)))


def f2():
    """ a, b, c三个点为等边三角形的三个顶点，它们距三角形中点为s,
    a 以速度v向b运动, b以速度v向c运动, c以速度v向a运动，问多长时间
    后三者相遇？
    """
    s = 1
    v = 1
    delta_t = 0.001
    t_total = 0
    a = [s, 0]
    b = [s, 120/180*np.pi]
    c = [s, 240/180*np.pi]
    fig = plt.subplot(polar=True)
    while s >0:
        plt.scatter([a[1], b[1], c[1]],
                    [a[0], b[0],c[0]],
                    linewidths=0)
        s = s-v*delta_t*np.cos(30/180*np.pi)
        a = [s, a[1]+v*delta_t*np.sin(30/180*np.pi)]
        b = [s, b[1]+v*delta_t*np.sin(30/180*np.pi)]
        c = [s, c[1]+v*delta_t*np.sin(30/180*np.pi)]
        t_total = t_total+delta_t
    plt.show()
    return t_total


if __name__ == '__main__':
    f1()
