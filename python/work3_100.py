#-------------------------------------------------------------------------------
# For the 3rd project
# By Dongjie, USTC/UM, start on Tue Sep 20 02:39:02 CST 2016
#-------------------------------------------------------------------------------

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
    bdate, edate = '2001-1-1', '2001-12-31'
    if True:
        omni = get_omni(bdate,edate, ['Bym', 'Bzm', 'V'], '1m')
        omni['EF'] = omni.V*(
                np.sqrt(omni.Bym*omni.Bym + omni.Bzm*omni.Bzm)-omni.Bzm)/2/1000
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

def func3():
    '''
    Case analysis: time constant of density response to IMF By
    '''
    #rtime = '2010-03-22 12:00:00' # quiet condition
    rtime = '2009-05-12 12:00:00' # quiet condition
    #rtime = '2010-10-17 06:00:00'
    rtime = pd.Timestamp(rtime)
    btime = rtime - pd.Timedelta('3 h')
    etime = rtime + pd.Timedelta('3 h')
    # IMF and Kan-Lee electric field variations
    def imf_ae():
        from matplotlib.ticker import AutoMinorLocator
        # By, Bz, K-L electric field
        fig1,ax1 = plt.subplots(3,1,sharex=True,figsize=(7,7))
        omni = get_omni(btime, etime, ['Bym', 'Bzm', 'V'], '1m')
        omni['EF'] = omni.V*(
                np.sqrt(omni.Bym*omni.Bym + omni.Bzm*omni.Bzm)-omni.Bzm)/2/1000
        pv = ['Bym', 'Bzm', 'EF']
        ylb = ['$B_y$ (nT)', '$B_z$ (nT)', 'Electric Field (mV/m)']
        yl = [(-20, 20), (-20, 20), (0, 20)]
        yt = [np.arange(-100, 101, 5),
              np.arange(-100, 101, 5), np.arange(0, 100, 2)]
        for k0 in range(3):
            plt.sca(ax1[k0])
            tmp = omni[pv[k0]]
            plt.plot((tmp.index-rtime)/pd.Timedelta('1h'), tmp, 'k')
            plt.xticks(np.arange(-30, 31, 1))
            plt.xlim((btime-rtime)/pd.Timedelta('1h'),
                     (etime-rtime)/pd.Timedelta('1h'))
            plt.ylabel(ylb[k0])
            plt.yticks(yt[k0])
            plt.ylim(yl[k0])
            plt.grid(dashes=(4,1))
            if k0 != 2:
                plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
        ax1[2].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
    imf_ae()

    # Read satellite data
    from scipy.interpolate import NearestNDInterpolator
    # arglats are np.nan near data head and tail, so the time interval of
    # data is `dt` longer than what is expected
    dt = pd.Timedelta('5h')
    denchamp = ChampDensity(btime-dt, etime+dt, satellite='champ')
    dengrace = ChampDensity(btime-dt, etime+dt, satellite='grace')
    dengoce = GoceData(btime-dt, etime+dt)
    if not dengoce.empty:
        from apexpy import Apex as apex
        mc = apex()
        dengoce['Mlat'], dengoce['MLT'] = mc.convert(
                dengoce['lat'], dengoce['long'], source='geo',
                dest='mlt', datetime=dengoce.index)
    den = [denchamp, dengrace, dengoce]
    fig2, ax2 = plt.subplots(2, 1, sharex=True)  # Test arglat and orbitn
    fig6, ax6 = plt.subplots(2, 1, sharex=True)  # Test altitude
    title = ('CHAMP', 'GRACE')#, 'GOCE')
    for k0 in range(2): # 0: Champ, 1: Grace, 2: Goce
        dent = den[k0]
        if dent.empty:
            continue
        dent['arglat'] = mf.lat2arglat(dent['lat'])
        dent['epoch'] = (dent.index - rtime)/pd.Timedelta('1h')
        # Calculate orbit number. It's better if no data gaps exist.
        orbitn = np.where(dent['arglat'].diff()<0, 1, 0)
        orbitn = orbitn.cumsum()
        dent['orbitn'] = orbitn

        # Test arglat and orbitn
        ax2[k0].plot(dent.epoch, dent.orbitn, 'k')
        ax2[k0].set_ylabel('Orbit Number (k)')
        ax21 = ax2[k0].twinx()
        ax21.plot(dent.epoch, dent.arglat, 'go')
        ax21.plot(dent.epoch, dent.lat, 'yo', alpha=0.3)
        ax21.set_xticks(np.arange(-30, 31, 1))
        ax21.set_xlim((btime-rtime)/pd.Timedelta('1h'),
                      (etime-rtime)/pd.Timedelta('1h'))
        ax21.set_ylabel('Arglat (g), Lat(y)')
        ax2[k0].set_title(title[k0])

        # Test altitude
        height = 'height' if k0<2 else 'alt'
        ax6[k0].plot(dent.epoch, dent[height], 'k')
        ax6[k0].set_ylabel('Height (km)')
        ax6[k0].set_xticks(np.arange(-30, 31, 1))
        ax6[k0].set_xlim((btime-rtime)/pd.Timedelta('1h'),
                      (etime-rtime)/pd.Timedelta('1h'))
        ax6[k0].set_title(title[k0])

        # Calculate delta density: rho(pass_n)-rho(pass_n-1)
        dent['ptime'] = np.nan
        dent['ddmin'] = np.nan
        for k1 in range(1,dent['orbitn'].max()+1):
            fpp = (dent['orbitn'] == k1-1)
            fpc = (dent['orbitn'] == k1)
            x0lat = dent.loc[fpp, 'lat'].copy()
            x0lon = dent.loc[fpp, 'LT'].copy()/12*180
            x1lat = dent.loc[fpc, 'lat'].copy()
            x1lon = dent.loc[fpc, 'LT'].copy()/12*180
            print(k0, x0lat.shape, x1lat.shape)
            dd = np.ones([x1lat.shape[0], x0lat.shape[0]])*np.nan
            for k2 in range(x1lat.shape[0]):
                dd[k2,:] = mf.great_circle_distance_earth(
                        x0lat, x0lon, x1lat.iloc[k2], x1lon.iloc[k2])
            ddminindex = np.argmin(dd, axis=1)
            dent.loc[fpc, 'ddmin'] = dd[range(dd.shape[0]), ddminindex]
            dent.loc[fpc, 'ptime'] = dent.index[fpp][ddminindex]
        pmsis = np.array(dent.loc[dent['ptime'], 'msis_rho'])
        cmsis = np.array(dent['msis_rho'])
        dent['pmsis_rho'] = pmsis
        dent['dmsis_rho'] = 100*(cmsis-pmsis)/pmsis
        prho400 = np.array(dent.loc[dent['ptime'], 'rho400'])
        crho400 = np.array(dent['rho400'])
        dent['prho400'] = prho400
        dent['drho400'] = 100*(crho400-prho400)/prho400
        prho = np.array(dent.loc[dent['ptime'], 'rho'])
        crho = np.array(dent['rho'])
        dent['prho'] = prho
        dent['drho'] = 100*(crho-prho)/prho
        #den[k0] = den[k0][btime:etime]
    ax2[-1].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
    ax6[-1].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))

    # mlat-mlt distribution of satellite orbits
    fig3, ax3 = plt.subplots(
            2,1,subplot_kw={'projection':'polar'}, figsize=[4.4,7])
    cl = list('rbk')
    lb = ['CHAMP', 'GRACE', 'GOCE']
    cls = [ChampDensity, ChampDensity, GoceData]
    for k11, k1 in enumerate(['N', 'S']):
        fc = 1 if k1 is 'N' else -1
        plt.sca(ax3[k11])
        for k0 in range(2): # Only CHAMP and GRACE
            dent = den[k0][btime:etime]
            if dent.empty:
                continue
            dent.__class__ = cls[k0]
            hc = dent.satellite_position_lt_lat(mag=True, ns=k1)
            hc.set_color(cl[k0])
            hc.set_label(lb[k0])
        mf.set_lat_lt_polar(plt.gca(), ns=k1, boundinglat=fc*0,
                            latgrids = np.arange(-60, 90, 30))
    hl = ax3[0].legend(loc=(0.85,0.8),fontsize=12)

    # Test ddmin
    fig4, ax4 = plt.subplots(2,1,sharex=True)
    for k0 in range(2):
        dent = den[k0].copy()
        if dent.empty:
            continue
        ax4[k0].plot(dent['epoch'], dent.ddmin, label=lb[k0], color=cl[k0])
        ax4[k0].set_xticks(np.arange(-30, 31, 1))
        ax4[k0].set_xlim((btime-rtime)/pd.Timedelta('1h'),
                         (etime-rtime)/pd.Timedelta('1h'))
        ax4[k0].set_title(lb[k0])
        ax4[k0].set_ylabel('$\Delta$r (km)')
    ax4[1].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))

    #  Density difference from orbit to orbit
    fig5, ax5 = plt.subplots(4,2,sharex=True, figsize=(8,9))
    title = ('CHAMP', 'GRACE')
    yl = (r'$\rho$ 400 km $(kg/m^{-3})$', r'$\Delta\rho$ 400 km (%)',
          r'Msis $\rho (kg/m^{-3})$', r'Msis $\Delta\rho$ (%)')
    col = ('rho400', 'drho400', 'msis_rho', 'dmsis_rho')
    for k0 in range(2):
        dent = den[k0].copy()
        if dent.empty:
            continue
        for k1 in range(4):
            ax5[k1, k0].plot(
                    dent.epoch, dent[col[k1]], label=lb[k0], color=cl[k0])
            dentmptmp = dent[dent.Mlat>0]
            ax5[k1, k0].plot(
                    dentmptmp.epoch, dentmptmp[col[k1]], 'ko', markersize=3)
            ax5[k1, k0].set_xticks(np.arange(-30, 31, 1))
            ax5[k1, k0].set_xlim((btime-rtime)/pd.Timedelta('1h'),
                             (etime-rtime)/pd.Timedelta('1h'))
            ax5[0, k0].set_title(title[k0])
            ax5[k1, 0].set_ylabel(yl[k1])
            if k1%2 ==1:
                ax5[k1, k0].axhline(color='gray', zorder=-1)
            ax5[-1,k0].set_xlabel(
                    'Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
    plt.tight_layout()

    # Density change as a function of argument of latitude
    fig7, ax7 = plt.subplots(3,2,sharex=True, figsize=(8,9))
    title = ('CHAMP', 'GRACE')
    yl = (r'$\rho$ $(kg/m^{-3})$',
          r'$\rho$ 400 km $(kg/m^{-3})$', r'Msis $\rho (kg/m^{-3})$')
    col = ('rho', 'rho400', 'msis_rho')
    for k0 in range(2):
        dent = den[k0][btime:etime]
        if dent.empty:
            continue
        for k1 in range(3):
            for k2 in range(dent['orbitn'].max()+1):
                dentt = dent[dent.orbitn==k2]
                ax7[k1, k0].plot(
                        dentt.arglat, dentt[col[k1]],
                        label=lb[k0], color=cl[k0])
                denttt = dentt[dentt.Mlat>0]
                ax7[k1, k0].plot(
                        denttt.arglat, denttt[col[k1]], 'ko', markersize=3)
            ax7[k1, k0].set_xticks(np.arange(0, 361, 60))
            ax7[k1, k0].set_xlim(0, 360)
            ax7[0, k0].set_title(title[k0])
            ax7[k1, 0].set_ylabel(yl[k1])
            ax7[-1,k0].set_xlabel('Argument of Latitude (degree)')
    plt.tight_layout()

    # Density change of the current and previous orbits
    fig8, ax8 = plt.subplots(3,2,sharex=True, figsize=(8,9))
    title = ('CHAMP', 'GRACE')
    yl = (r'$\rho$ $(kg/m^{-3})$',
          r'$\rho$ 400 km $(kg/m^{-3})$', r'Msis $\rho (kg/m^{-3})$')
    col = (('rho', 'prho'), ('rho400', 'prho400'), ('msis_rho', 'pmsis_rho'))
    for k0 in range(2):
        dent = den[k0][btime:etime]
        if dent.empty:
            continue
        for k1 in range(3):
            ax8[k1, k0].plot(
                    dent.epoch, dent[col[k1][0]],
                    label=lb[k0], color=cl[k0])
            ax8[k1, k0].plot(
                    dent.epoch, dent[col[k1][1]],
                    label=lb[k0], color='k', zorder=0)
            ax8[k1, k0].set_xticks(np.arange(-30, 31, 1))
            ax8[k1, k0].set_xlim((btime-rtime)/pd.Timedelta('1h'),
                             (etime-rtime)/pd.Timedelta('1h'))
            ax8[0, k0].set_title(title[k0])
            ax8[k1, 0].set_ylabel(yl[k1])
            ax8[-1,k0].set_xlabel(
                    'Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
    plt.tight_layout()
    #    fig1.savefig('/home/guod/Documents/work/fig/w03_func03_02a')
    #    fig2.savefig('/home/guod/Documents/work/fig/w03_func03_02b')
    #    fig3.savefig('/home/guod/Documents/work/fig/w03_func03_02c')
    #    fig4.savefig('/home/guod/Documents/work/fig/w03_func03_02d')
    #    fig5.savefig('/home/guod/Documents/work/fig/w03_func03_02e')
    #    fig6.savefig('/home/guod/Documents/work/fig/w03_func03_02f')
    #    fig7.savefig('/home/guod/Documents/work/fig/w03_func03_02g')
    #    fig8.savefig('/home/guod/Documents/work/fig/w03_func03_02h')
    return
# END
#------------------------------------------------------------
if __name__ == '__main__':
    plt.close('all')
    a = func3()
    plt.show()
