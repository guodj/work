#!/home/gdj/anaconda3/bin/python3
#-*- coding: utf-8 -*-

__author__ = 'Guo Dongjie'

"""
Some funny and useful python programmes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

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


def f4():
    """买房与租房
    """
    area = 100
    # 买房
    bp = 10000  #每平米房价（没算装修，税等）
    mp = 100000 #装修费用
    r3 = 0.3   #首付比例
    r4 = 0.049 #商贷利率
    r5 = 0.0325 #公积金贷款利率
    p1 = 550000    #公积金贷款限额
    y = 20    #贷款年限
    r1 = 0.02  #预期房价每年上升率

    # 租房
    rp = 2500  #假设为同等类型的住房的每月租价
    r2 = 0.01  #预期房租每年上升率
    r6 = 0.05  #投资预期年收益

    shoufu = bp*area*r3  #首付总数
    if bp*area*(1-r3)<=p1:
        fangdai1 = bp*area*(1-r3)*(r5/12+(r5/12)/((r5/12+1)**(12*y)-1))
        fangdai2 = 0
    else:
        fangdai1 = p1*(r5/12+(r5/12)/((r5/12+1)**(12*y)-1))
        fangdai2 = (bp*area*(1-r3)-p1)*(r4/12+(r4/12)/((r4/12+1)**(12*y)-1))
    fangdai = fangdai1+fangdai2
    if fangdai<rp:
        print('房贷低于房租，肯定买房好')
        return
    mfztr = fangdai*12*y+shoufu  #买房总投入
    fwjz = bp*area*(1+r1)**y   #y年后预期房价
    tzzsr = shoufu*(1+r6)**y + (fangdai-rp)*12*((1+r6)**y-1)/r6 #如果租房可以获得的投资总收入

    print('租房VS买房\n--------------------------------------------------')
    print('假设:\n当前房价为{:.2f}万/m^2，并且预期每年上涨{:.2f}%.'.format(bp/1e4, r1*100))
    print('当前政策为首付比例{:.0f}%, 商贷利率{:.2f}%，'
          '公积金贷款利率为{:.2f}%，公积金贷款限额{:.2f}万。'.
          format(r3*100,r4*100,r5*100,p1/1e4))
    print('--------------------------------------------------')
    print('若买房（{:d}m^2，贷款{:d}年）：'.format(area,y))
    print('每月房贷： {:.2f}元'.format(fangdai))
    print('共计投入资金（首付+贷款）：{:.2f}万'.format(mfztr/1e4))
    print('买房总收入（房屋{:d}年后价值-总投入）： {:.2f}万。'.format(y,(fwjz-mfztr)/1e4))
    print('--------------------------------------------------')
    print('若租房:')
    print('设租房每月租金{:d}元，首付与剩余租金用于投资，预期年收益为{:.2f}%.'.
          format(rp,r6*100))
    print('租房总收入：{:.2f}万。'.format((tzzsr-mfztr)/1e4))
    print('--------------------------------------------------')
    print('注意：')
    print('未考虑租金上涨；未考虑装修；'
          '未考虑中国租房制度对租房质量的影响；未考虑房子与户口，医疗，孩子教育的关系。')

def f5(x=1000,r=0.1):
    """每月投资x,投资n年，可以有多少钱？
    """
    y1 = np.sum(np.arange(1,13)/12)
    print('每月投资{:d}元,设年收益{:.0f}%'.format(x,r*100))
    for n in np.arange(5,31,5):
        print('{:d}年后，将收获{:.2f}万'.format(n,x*(12+r*y1)*((1+r)**n-1)/r/1e4))
    return


def f6():
    """定投沪深300指数：
    注：不考虑管理费，红利等
    注：定投策略的前提是指数长期见长
    """
    myfont = matplotlib.font_manager.FontProperties(fname='/home/gdj/.local/share/fonts/simsun.ttc')
    #myfont = matplotlib.font_manager.FontProperties(fname='/home/gdj/.fonts/s/SIMSUN.ttc')
    matplotlib.rcParams['axes.unicode_minus'] = False
    hs300 = pd.read_csv(
            '/data/hs300.csv',encoding='gbk',
            skiprows=1,
            header=None,
            names=['date','closep','maxp','minp','openp','volume','amo'],
            usecols=[0,3,4,5,6,7,8],
            index_col=0,
            parse_dates=0)
    hs300 = hs300.sort_index()
    hs300m = pd.DataFrame()
    hs300m['openp'] = hs300['openp'].resample('1M').first()
    hs300m['closep'] = hs300['closep'].resample('1M').last()
    hs300m['maxp'] = hs300['maxp'].resample('1M').max()
    hs300m['minp'] = hs300['minp'].resample('1M').min()
    hs300m['volume'] = hs300['volume'].resample('1M').sum()
    hs300m['amo'] = hs300['amo'].resample('1M').sum()
    #hs300m = hs300m['2010-1'::]

    plt.figure()
    # 币值成本平均定投
    ii = 1000 # 初始投资金额
    r = 0.008 # 每月增加的投资金额
    t = np.arange(len(hs300m)) #总月数为len(hs300m)
    cb = ii*(1+r)**t #每月成本
    zcb = np.cumsum(cb) # 总成本
    v = zcb*np.nan
    v[0] = ii
    for k in np.arange(1,len(hs300m)):
        pt = hs300m.ix[k,'closep']
        pl = hs300m.ix[k-1,'closep']
        v[k] = v[k-1]*pt/pl+cb[k]
    plt.plot(hs300m.index,v,'b--',label='收入：币值成本平均')
    #plt.plot(hs300m.index,zcb,'r--',label='成本：币值成本平均')
    #plt.plot(hs300m.index,v-zcb,'b-',label='净收入：币值成本平均')
    plt.grid()

    #价值平均定投
    ii = 1000  #第一次投资
    r = 0.008 #月增长率
    t = np.arange(len(hs300m)) #总月数为len(hs300m)
    v = ii*(t+1)*(1+r)**t  #价值曲线
    cb = v*np.nan  #每月成本,指的是单纯工资支出
    cb[0] = ii
    zcb = v*np.nan  #总成本,指的是单纯工资总支出
    zcb[0] = ii
    tq = v*np.nan  #由于基金净值增加，从基金中提取的钱存到此账户(账户余额>0)
    tq[0] = 0
    for k in np.arange(1,len(cb)):
        nm = v[k]-v[k-1]*hs300m.ix[k,'closep']/hs300m.ix[k-1,'closep']
        tq[k] = tq[k-1]-nm
        if tq[k]<0:
            cb[k] = -tq[k]
            tq[k] = 0
        else:
            cb[k] = 0
        zcb[k] = zcb[k-1]+cb[k]
    plt.plot(hs300m.index,v,'b',label='收入：价值平均')
    plt.plot(hs300m.index,tq,'r',label='收入：价值平均')
    plt.plot(hs300m.index,zcb,'k',label='成本：价值平均')
    #plt.plot(hs300m.index,v+tq-zcb,'r-',label='净收入：价值平均')
    plt.legend(prop=myfont,loc=2)
    plt.figure()
    plt.plot(hs300.index,hs300.closep,'gray',marker='o',markersize=3,label='沪深300指数')
    plt.show()
    return hs300


def f7():
    """沪深300指数每月的哪一天最小：
    """
    myfont = matplotlib.font_manager.FontProperties(fname='/home/gdj/.local/share/fonts/simsun.ttc')
    #myfont = matplotlib.font_manager.FontProperties(fname='/home/gdj/.fonts/s/SIMSUN.ttc')
    matplotlib.rcParams['axes.unicode_minus'] = False
    hs300 = pd.read_csv(
            '/data/hs300.csv',encoding='gbk',
            skiprows=1,
            header=None,
            names=['date','closep','maxp','minp','openp','volume','amo'],
            usecols=[0,3,4,5,6,7,8],
            index_col=0,
            parse_dates=0)
    hs300 = hs300.sort_index()
    hs300 = hs300[hs300.volume!=0]
    hs3001 = hs300.groupby([hs300.index.year,hs300.index.day])['closep','maxp','minp','openp'].median()
    fig = plt.figure()
    for k0 in hs3001.index.get_level_values(0):
        plt.plot(hs3001.loc[k0].index,hs3001.loc[k0].closep,'b',marker='o')
    plt.xlim(1,31)
    #没有统计上的规律
    plt.show()

    return hs3001


if __name__ == '__main__':
    plt.close('all')
    a=f7()
