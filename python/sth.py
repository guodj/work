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
    glon, glat = np.meshgrid(np.arange(-180,180,1),np.arange(-90,91,1))
    mlat, mlon = m.convert(glat, glon, 'geo', 'apex', height=0,datetime=pd.Timestamp('2016-5-7'))
    mlat1, mlt = m.convert(glat, glon, 'geo', 'mlt', datetime=pd.Timestamp('2016-5-7'))

    fig = plt.figure(figsize=(7,11))

    ax = plt.subplot(4,2,1)
    mp = Basemap(projection='npstere',boundinglat=20,lon_0=0,resolution='l',round=True)
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-80,81,20.),dashes=(10,1),linewidth=2)
    mp.drawmeridians(np.arange(0,360,60.),labels=[True,True,True,True],
                     dashes=(10,1),linewidth=2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlon,latlon=True,levels=np.arange(-180,180,30),colors='r')
    plt.clabel(cs1,np.arange(0,81,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('NH',loc='left',color='r')

    ax = plt.subplot(4,2,2)
    mp = Basemap(projection='spstere',boundinglat=-20,lon_0=0,resolution='l',round=True)
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-80,81,20.),dashes=(10,1),linewidth=2)
    mp.drawmeridians(np.arange(0,360,60.),labels=[True,True,True,True],
                     dashes=(10,1),linewidth=2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlon,latlon=True,levels=np.arange(-180,180,30),colors='r')
    plt.clabel(cs1,np.arange(-80,1,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('SH',loc='left',color='r')

    ax = plt.subplot(4,1,2)
    mp = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-90,
                 urcrnrlon=180,urcrnrlat=90,resolution='l')
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-90,91,30.),dashes=(10,1),linewidth=2,
                     labels=[True,False,False,False])
    mp.drawmeridians(np.arange(0,360,60.),labels=[False,False,False,True],
                     dashes=(10,1),linewidth=2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlon,latlon=True,levels=np.arange(-180,180,30),colors='r')
    plt.clabel(cs1,np.arange(-80,81,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('apex in geodetic')

    ax = plt.subplot(4,1,3)
    mp = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-90,
                 urcrnrlon=180,urcrnrlat=90,resolution='l')
    mp.drawcoastlines(linestyle='--',color='black')
    mp.drawparallels(np.arange(-90,91,30.),dashes=(10,1),linewidth=2,
                     labels=[True,False,False,False])
    mp.drawmeridians(np.arange(0,360,60.),labels=[False,False,False,True],
                     dashes=(10,1),linewidth=2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(glon,glat,mlat1,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(glon,glat,mlt,latlon=True,levels=np.arange(0,24,2),colors='r')
    plt.clabel(cs1,np.arange(-80,81,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('mlat, mlt in geodetic')

    mlon, mlat = np.meshgrid(np.arange(-180,180,1),np.arange(-90,91,1))
    glat, glon = m.convert(mlat, mlon, 'apex', 'geo', height=0,datetime=pd.Timestamp('2016-5-7'))
    ax = plt.subplot(4,1,4)
    mp = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-90,
                 urcrnrlon=180,urcrnrlat=90,resolution='l')
    mp.drawparallels(np.arange(-90,91,30.),dashes=(10,1),linewidth=2,
                     labels=[True,False,False,False])
    mp.drawmeridians(np.arange(0,360,60.),labels=[False,False,False,True],
                     dashes=(10,1),linewidth=2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cs1 = mp.contour(mlon,mlat,glat,latlon=True,levels=np.arange(-80,81,20),colors='r')
    cs2 = mp.contour(mlon,mlat,glon,latlon=True,levels=np.arange(-180,181,30),colors='r')
    plt.clabel(cs1,np.arange(-80,81,20),fontsize=12, inline=True,fmt='%d',colors='b')
    plt.title('geodetic in apex')
    plt.subplots_adjust(left=0.1,bottom=0.05,hspace=0.3)
    plt.show()
    return mlat, mlon


if __name__ == '__main__':
    a=f3()
