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


def f3():
    """地球磁场的位形
    """
    from apexpy import Apex
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.basemap import Basemap

    m = Apex()
    glon, glat = np.meshgrid(np.arange(-180,180.1,5),np.arange(-90,90.1,5))
    mlat, mlon = m.convert(glat, glon, 'geo', 'apex', height=0)
    mlat1, mlt = m.convert(glat, glon, 'geo', 'mlt', datetime=pd.Timestamp('2016-5-7'))

    fig, ax = plt.subplots(1,2,figsize=(13.5,6))
    plt.sca(ax[0])
    mp = Basemap(projection='npstere',boundinglat=1,lon_0=0,resolution='l',round=True)
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-80.,81.,20.),labels=[False,True,True,False])
    mp.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,True,True])
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlon,latlon=True,levels=np.arange(-180,181,20),colors='r')
    plt.clabel(cs1,np.arange(0,81,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('NH')

    plt.sca(ax[1])
    mp = Basemap(projection='spstere',boundinglat=-1,lon_0=0,resolution='l',round=True)
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-80.,81.,20.))
    mp.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,True,True])
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlon,latlon=True,levels=np.arange(-180,181,20),colors='r')
    plt.clabel(cs1,np.arange(-80,1,20),fontsize=9, inline=True,fmt='%d',colors='b')
    plt.title('SH')
    plt.show()
    return mlat, mlon


if __name__ == '__main__':
    a=f3()
