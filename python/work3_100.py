#--------------------------------------------------------------------------------
# For the 3rd project
#
# By Dongjie, USTC/UM, start on Tue Sep 20 02:39:02 CST 2016
#--------------------------------------------------------------------------------

#Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from champ_grace import *
from goce import *
from omni import *
import myfunctions as mf

DATADIR = '/home/guod/data/'
#DATADIR = '/data/'
def func1():
    '''
    Find interesting cases in 2009 (11,12) and 2010
    interesting case:
        1, Sharp IMF By reversion
        2, Small Bz and electric field change
    '''
    # For IMF and Kan-Lee electric field
    from pylab import draw # Use draw()
    import matplotlib.dates as mdates
    from matplotlib.ticker import AutoMinorLocator
    bdate, edate = '2002-1-1', '2002-12-31'
    if False:
        omni = get_omni(bdate,edate, ['Bym', 'Bzm', 'V'], '1m')
        omni['EF'] = omni.V*(np.sqrt(omni.Bym*omni.Bym + omni.Bzm*omni.Bzm)-omni.Bzm)/2/1000
        pd.to_pickle(omni,DATADIR + 'tmp/w3_01_1.dat')
    omni = pd.read_pickle(DATADIR + 'tmp/w3_01_1.dat')
    fig1,ax1 = plt.subplots(3,1,sharex=True,figsize=(7,7)) # By, Bz, E
    pv = ['Bym', 'Bzm', 'EF']
    ylb = ['$B_y$ (nT)', '$B_z$ (nT)', 'Electric Field (mV/m)']
    yl = [(-20, 20), (-20, 20), (0, 5)]
    yt = [np.arange(-20, 21, 5), np.arange(-20, 21, 5), np.arange(0, 6, 1)]
    hours = mdates.HourLocator(range(0,25,3))
    hoursfmt = mdates.DateFormatter('%H')
    for k0 in range(3):
        plt.sca(ax1[k0])
        tmp = omni[pv[k0]][bdate:edate]
        plt.plot(tmp.index, tmp)
        plt.xlim(bdate, edate)
        plt.ylabel(ylb[k0])
        plt.yticks(yt[k0])
        plt.ylim(yl[k0])
        plt.grid(dashes=(4,1))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
    for k0 in range(2):
        plt.sca(ax1[k0])
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
    plt.sca(ax1[2])
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(1))
    plt.xlabel('Hours of {:s}'.format(pd.Timestamp(bdate).strftime('%Y-%m-%d')))
    for date in pd.date_range(bdate,edate):
        ax1[-1].set_xlim(date,date+pd.Timedelta('2D'))
        ax1[-1].xaxis.set_major_locator(hours)
        ax1[-1].xaxis.set_major_formatter(hoursfmt)
        ax1[-1].set_xlabel('Hours of date: '+
                           date.date().strftime('%Y-%m-%d')+
                           '/'+(date+pd.Timedelta('1D')).date().strftime('%d'))
        # By default, figures won't change until end of the script
        # draw() forces a figure redraw
        draw()
        plt.show()
        input()
    return

def func2():
    #----------------------------------------
    # Case 1
    #    bdate = '2010-5-29 03:00:00'
    #    mdate = '2010-5-29 12:00:00'
    #    edate = '2010-5-29 21:00:00'
    # Case 2
    #    bdate = '2009-11-14 06:00:00'
    #    mdate = '2009-11-14 11:00:00'
    #    edate = '2009-11-14 18:00:00'
    # Case 3
    #    bdate = '2009-12-12 15:00:00'
    #    mdate = '2009-12-12 20:00:00'
    #    edate = '2009-12-13 06:00:00'
    # Case 4
    #    bdate = '2010-01-02 02:00:00'
    #    mdate = '2010-01-02 07:00:00'
    #    edate = '2010-01-02 12:00:00'
    # Case 5
    #    bdate = '2010-01-03 03:00:00'
    #    mdate = '2010-01-03 12:00:00'
    #    edate = '2010-01-03 16:00:00'
    # Case 6
    #    bdate = '2010-01-05 09:00:00'
    #    mdate = '2010-01-05 18:00:00'
    #    edate = '2010-01-06 03:00:00'
    # Case 7
    #    bdate = '2010-04-12 09:00:00'
    #    mdate = '2010-04-12 18:00:00'
    #    edate = '2010-04-13 00:00:00'
    # Case 8
    bdate = '2010-05-18 21:00:00'
    mdate = '2010-05-19 06:00:00'
    edate = '2010-05-19 23:00:00'
    #----------------------------------------
    dench = ChampDensity(bdate,edate,satellite='champ')
    dench['arglat'] = mf.lat2arglat(dench.lat)
    mdench = dench.groupby([np.floor(dench.arglat/3)*3, dench.index<mdate])['rho'].mean()
    dengr = ChampDensity(bdate,edate,satellite='grace')
    dengr['arglat'] = mf.lat2arglat(dengr.lat)
    mdengr = dengr.groupby([np.floor(dengr.arglat/3)*3, dengr.index<mdate])['rho'].mean()
    dengo = GoceData(bdate,edate)
    dengo['arglat'] = mf.lat2arglat(dengo.lat)
    mdengo = dengo.groupby([np.floor(dengo.arglat/3)*3, dengo.index<mdate])['rho'].mean()
    den = (dengo, dench, dengr)
    mden = (mdengo, mdench, mdengr)
    fig, ax = plt.subplots(3,2,sharex='col', sharey='row', figsize=(7.4,7.9))
    blat = 60
    xl = ((blat, 180-blat), (180+blat, 360-blat))
    yl = ((5, 35), (0, 20), (0, 0.3))
    ylb = ('GOCE, $10^{-12} kg/m^{-3}$',
           'CHAMP, $10^{-12} kg/m^{-3}$',
           'GRACE, $10^{-12} kg/m^{-3}$')
    tl = ('North', 'South')
    for k0 in range(3):
        for k1 in range(2):
            plt.sca(ax[k0,k1])
            tmp1, tmp2 = den[k0], mden[k0]
            plt.plot(tmp1[bdate:mdate].arglat, tmp1[bdate:mdate].rho/1e-12,
                     'o', color='lightcoral', alpha=1, markersize=3)
            plt.plot(tmp1[mdate:edate].arglat, tmp1[mdate:edate].rho/1e-12,
                     'o', color='lightblue', alpha=1, markersize=3)
            line1, = plt.plot(tmp2[:,False].index, tmp2[:,False]/1e-12, 'b', label='After')
            line2, = plt.plot(tmp2[:,True].index, tmp2[:,True]/1e-12, 'r', label='Before')
            plt.xticks(np.arange(0,361,10))
            plt.xlim(xl[k1])
            plt.ylim(yl[k0])
            plt.gca().set_frame_on(True)
            if k0 == 0:
                plt.title(tl[k1])
            if k1 == 0:
                plt.ylabel(ylb[k0])
            if (k0 == 2) & (k1 == 1):
                plt.legend(handles=(line2, line1))
            if k0 == 2:
                plt.xlabel('Argument of Latitude')
            plt.grid('on')
    plt.tight_layout()
    return

def func3():
    '''
    Case analysis: time constant of density response to IMF By
    '''
    rtime = '2009-07-05 20:50:00'
    rtime = pd.Timestamp(rtime)
    btime = rtime - pd.Timedelta('1.5 h')
    etime = rtime + pd.Timedelta('3 h')
    # IMF and Kan-Lee electric field variations
    from matplotlib.ticker import AutoMinorLocator
    fig1,ax1 = plt.subplots(3,1,sharex=True,figsize=(7,7)) # By, Bz, Kan-Lee electric field
    omni = get_omni(btime,etime, ['Bym', 'Bzm', 'V'], '1m')
    omni['EF'] = omni.V*(np.sqrt(omni.Bym*omni.Bym + omni.Bzm*omni.Bzm)-omni.Bzm)/2/1000
    pv = ['Bym', 'Bzm', 'EF']
    ylb = ['$B_y$ (nT)', '$B_z$ (nT)', 'Electric Field (mV/m)']
    yl = [(-20, 20), (-20, 20), (0, 10)]
    yt = [np.arange(-100, 101, 5), np.arange(-100, 101, 5), np.arange(0, 100, 1)]
    for k0 in range(3):
        plt.sca(ax1[k0])
        tmp = omni[pv[k0]]
        plt.plot((tmp.index-rtime)/pd.Timedelta('1m'), tmp)
        plt.xlim((btime-rtime)/pd.Timedelta('1m'), (etime-rtime)/pd.Timedelta('1m'))
        plt.ylabel(ylb[k0])
        plt.yticks(yt[k0])
        plt.ylim(yl[k0])
        plt.grid(dashes=(4,1))
        if k0 != 2:
            plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
    # mlat-mlt distribution of satellite orbits
    fig2, ax2 = plt.subplots(2,1,subplot_kw={'projection':'polar'}, figsize=[4.4,7])
    dt = pd.Timedelta('5h')
    denchamp = ChampDensity(btime-dt,etime+dt,satellite='champ')
    dengrace = ChampDensity(btime-dt,etime+dt,satellite='grace')
    dengoce = GoceData(btime,etime)
    den = [denchamp, dengrace, dengoce]
    for k0 in range(3):
        if den[k0].empty:
            continue
        den[k0]['arglat'] = mf.lat2arglat(den[k0]['lat'])
        den[k0]['epoch'] = (den[k0].index - rtime)/pd.Timedelta('1m')
        den[k0] = den[k0][btime:etime]
    cl = list('rbk')
    lb = ['CHAMP', 'GRACE', 'GOCE']
    cls = [ChampDensity, ChampDensity, GoceData]
    for k11, k1 in enumerate(['N', 'S']):
        fc = 1 if k1 is 'N' else -1
        plt.sca(ax2[k11])
        for k0 in range(3):
            dentmp = den[k0]
            if dentmp.empty:
                continue
            dentmp.__class__ = cls[k0]
            hc = dentmp.satellite_position_lt_lat(mag=True, ns=k1)
            hc.set_color(cl[k0])
            hc.set_label(lb[k0])
        mf.set_lat_lt_polar(plt.gca(), ns=k1, boundinglat=fc*60)
    hl = ax2[0].legend(loc=(0.85,0.8),fontsize=12)
    # delta density from reversing time of IMF By
    from scipy.interpolate import interp1d
    fig = plt.figure()
    fden = [None, None, None]
    for k0 in range(3):
        dentmp = den[k0].copy()
        if dentmp.empty:
            continue
        arglat0 = dentmp.loc[btime:rtime, 'arglat']
        rho0 = dentmp.loc[btime:rtime, 'rho400']
        ff = interp1d( arglat0, rho0, kind='linear', fill_value='extrapolate')
        arglat1 = dentmp['arglat']
        dentmp['rhoitp'] = ff(arglat1)
        dentmp = dentmp[dentmp.lat<-60]
        plt.plot(dentmp['epoch'], 100*(dentmp.rho400-dentmp.rhoitp)/dentmp.rhoitp, label=lb[k0])
        plt.xlim(dentmp['epoch'].min(), dentmp['epoch'].max())
    plt.legend()
    return

# END
#------------------------------------------------------------
if __name__ == '__main__':
    plt.close('all')
    a = func3()
    plt.show()
