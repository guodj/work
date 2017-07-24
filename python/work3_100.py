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
import seaborn as sns
from apexpy import Apex
import os
import gitm
import gitm_3D_const_alt as g3ca
import matplotlib.animation as animation
import gitm_vorticity as gv
from matplotlib.ticker import AutoMinorLocator
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs
from pylab import draw # Use draw()
import matplotlib.dates as mdates
from scipy.interpolate import NearestNDInterpolator
from spacepy.datamodel import dmarray
import gitm_divergence as gd
from mpl_toolkits.axes_grid1 import AxesGrid
import gitm_gradient as ggt
sns.set('paper', 'whitegrid')

DATADIR = os.environ.get('DATAPATH')
Re = 6371*1000 # Earth radius, unit: m
def func1():
    '''
    Find interesting cases in 2009 (11,12) and 2010
    interesting case:
        1, Sharp IMF By reversion
        2, Small Bz and electric field change
    '''
    # For IMF and Kan-Lee electric field
    bdate, edate = '2010-1-1', '2010-12-31'
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
    #rtime = pd.Timestamp('2010-02-27 00:00:00') # quiet condition
    #rtime = pd.Timestamp('2009-05-12 12:00:00') # quiet condition
    rtime = pd.Timestamp('2004-07-25 07:00:00') # active condition
    #rtime = '2010-10-17 06:00:00'
    btime = rtime - pd.Timedelta('3 h')
    etime = rtime + pd.Timedelta('3 h')

    figdir = '/home/gdj/documents/work/fig/'
    if True:  # IMF and Kan-Lee electric field variations
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
        fig1.savefig(figdir + 'w03_func03_03a')

    lb = ['CHAMP', 'GRACE', 'GOCE']
    cl = list('rbk')
    if True:  # Read satellite data
        # arglats are np.nan near data head and tail, so the time interval of
        # data is `dt` longer than what is expected
        dt = pd.Timedelta('5h')
        denchamp = ChampDensity(btime-dt, etime+dt, satellite='champ')
        dengrace = ChampDensity(btime-dt, etime+dt, satellite='grace')
        dengoce = GoceData(btime-dt, etime+dt)
        if not dengoce.empty:
            mc = apex()
            dengoce['Mlat'], dengoce['MLT'] = mc.convert(
                    dengoce['lat'], dengoce['long'], source='geo',
                    dest='mlt', datetime=dengoce.index)
        den = [denchamp, dengrace, dengoce]
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
            dent['dmsis_rho'] = 100*(cmsis-pmsis)/pmsis  # Relative difference
            #dent['dmsis_rho'] = cmsis-pmsis # Absolute difference
            prho400 = np.array(dent.loc[dent['ptime'], 'rho400'])
            crho400 = np.array(dent['rho400'])
            dent['prho400'] = prho400
            # Relative difference
            dent['drho400'] = 100*(crho400-prho400)/prho400
            #dent['drho400'] = crho400-prho400  # Absolute difference
            prho = np.array(dent.loc[dent['ptime'], 'rho'])
            crho = np.array(dent['rho'])
            dent['prho'] = prho
            dent['drho'] = 100*(crho-prho)/prho  # Relative difference
            #dent['drho'] = crho-prho  # Absolute difference
            #den[k0] = den[k0][btime:etime]

    if True:  # Test arglat, orbitn and altitude
        fig2, ax2 = plt.subplots(2, 1, sharex=True)  # Test arglat and orbitn
        fig6, ax6 = plt.subplots(2, 1, sharex=True)  # Test altitude
        title = ('CHAMP', 'GRACE')#, 'GOCE')
        for k0 in range(2): # 0: Champ, 1: Grace, 2: Goce
            dent = den[k0]
            if dent.empty:
                continue
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
        ax2[-1].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
        ax6[-1].set_xlabel('Hours from '+rtime.strftime('%Y-%m-%d %H%M UT'))
        fig2.savefig(figdir + 'w03_func03_03b')
        fig6.savefig(figdir + 'w03_func03_03f')

    if True:  # mlat-mlt distribution of satellite orbits
        fig3, ax3 = plt.subplots(
                2,1,subplot_kw={'projection':'polar'}, figsize=[4.4,7])
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
        fig3.savefig(figdir + 'w03_func03_03c')

    if True:  # Test ddmin
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
        fig4.savefig(figdir + 'w03_func03_03d')

    if True:  #  Density difference from orbit to orbit
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
                dentmptmp = dent[dent.lat>0]
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
        fig5.savefig(figdir + 'w03_func03_03e')

    if True:  # Density change as a function of argument of latitude
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
        fig7.savefig(figdir + 'w03_func03_03g')

    if True:  # Density change of the current and previous orbits
        fig8, ax8 = plt.subplots(3,2,sharex=True, figsize=(8,9))
        title = ('CHAMP', 'GRACE')
        yl = (r'$\rho$ $(kg/m^{-3})$',
              r'$\rho$ 400 km $(kg/m^{-3})$', r'Msis $\rho (kg/m^{-3})$')
        col = (('rho', 'prho'), ('rho400', 'prho400'),
               ('msis_rho', 'pmsis_rho'))
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
        fig8.savefig(figdir + 'w03_func03_03h')
    return


def animate_all():
    stime = pd.Timestamp('2010-03-23 00:00:00')
    etime = pd.Timestamp('2010-03-23 05:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    timeidxepoch = (timeidx - pd.Timestamp('2010-03-23 00:00:00'))\
            /pd.Timedelta('1hour')
    fn1 = ['/home/guod/big/raid4/guod/run_imfby/run1c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]
    fn2 = ['/home/guod/big/raid4/guod/run_imfby/run2c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]

    # save path
    path =  '/home/guod/Documents/work/fig/w03/animation/'

    # geomagnetic pole
    apex = Apex(date=2010)
    qlatn, qlonn = apex.convert(90, 0, source='apex', dest='geo', height=400)
    qlats, qlons = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    fig = plt.figure(figsize=[8,16])
    def animate_den_wind(i):
        # IMF By
        plt.subplot(4, 1, 1)
        plt.gca().clear()
        btimeby = pd.Timestamp('2010-03-20 00:00:00')
        etimeby = pd.Timestamp('2010-03-24 00:00:00')
        dtby = pd.Timedelta('1min')
        timeby = pd.date_range(btimeby, etimeby, freq=dtby)
        epochhour = (timeby - pd.Timestamp('2010-03-23 00:00:00'))\
                /pd.Timedelta('1hour')
        fnimf1 = '/home/guod/WD4T/gitm/run_imfby/run1c/imf1.dat'
        imf1 = pd.read_csv(
                fnimf1, delim_whitespace=True, comment='#',
                header=None,
                names=('year', 'month', 'day', 'hour', 'minute', 'second',
                       'ms', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
                usecols=['by'])
        imf1 = pd.DataFrame(np.array(imf1), index=timeby, columns=['By'])
        plt.plot(epochhour, imf1.By, 'b')
        fnimf2 = '/home/guod/WD4T/gitm/run_imfby/run2c/imf2.dat'
        imf2 = pd.read_csv(
                fnimf2, delim_whitespace=True, comment='#',
                header=None,
                names=('year', 'month', 'day', 'hour', 'minute', 'second',
                       'ms', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
                usecols=['by'])
        imf2 = pd.DataFrame(np.array(imf2), index=timeby, columns=['By'])
        plt.plot(epochhour, imf2.By, 'r')
        plt.xlim(-1, 5)
        plt.xticks(np.arange(-1, 5.1, 0.5))
        plt.xlabel('Epoch (hour)')
        plt.ylim([-10, 10])
        plt.ylabel(r'IMF $B_Y$')
        plt.axvline(timeidxepoch[i], color='k', linestyle='--')

        # read gitm data
        g1, g2 = [gitm.GitmBin(
                k[i], varlist=['Rho', 'Temperature', 'V!Dn!N (north)',
                               'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                               'V!Di!N (north)']) for k in [fn1, fn2]]

        # create axis
        ax = list(range(6))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 50] if ins==0 else [-50, -90]
            for irun in range(3):
                ax[ins+irun*2], projection[ins+irun*2] = gcc.create_map(
                        4, 2, 3+ins+irun*2, 'polar', nlat=nlat, slat=slat,
                        dlat=10, centrallon=g3ca.calculate_centrallon(
                            g1, 'polar',  useLT=True),
                        coastlines=False)

        # Density
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=400)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2, alt=400)
        hc = [ax[k].contourf(
                lon1, lat1, zdata1, transform=ccrs.PlateCarree(),
                levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
                extend='both') for k in [0, 1]]
        hc = [ax[k].contourf(
                lon2, lat2, zdata2, transform=ccrs.PlateCarree(),
                levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
                extend='both') for k in [2, 3]]
        diffzdata = 100*(zdata2-zdata1)/zdata1
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, transform=ccrs.PlateCarree(),
                levels=np.linspace(-30, 30, 21), cmap='seismic',
                extend='both') for k in [4, 5]]

        # wind
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=400)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=400)
        for iax in range(6):
            if iax == 0 or iax == 1:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind1.copy(), nwind1.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 2 or iax == 3:
                lon0, lat0, ewind, nwind = (
                        lon2.copy(), lat2.copy(), ewind2.copy(), nwind2.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 4 or iax == 5:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            hq = ax[iax].quiver(
                    lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                    color='k', regrid_shape=20)
            # ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
            # hc = plt.colorbar(hc, ticks=np.arange(3, 7)*1e-12)
            # hc.set_label(r'$\rho$ (kg/m$^3$)')
            ax[iax].scatter(
                    qlonn, qlatn, color='k', transform=ccrs.PlateCarree())
            ax[iax].scatter(
                    qlons, qlats, color='k', transform=ccrs.PlateCarree())
        return
    anim = animation.FuncAnimation(
            fig, animate_den_wind, interval=200, frames=len(fn1))
    anim.save(path+'den_wind.gif', writer='imagemagick')
    plt.close('all')

    #    fig = plt.figure(figsize=[8,16])
    #    def animate_temperature(i):
    #        # IMF By
    #        plt.subplot(4, 1, 1)
    #        plt.gca().clear()
    #        btimeby = pd.Timestamp('2010-03-20 00:00:00')
    #        etimeby = pd.Timestamp('2010-03-24 00:00:00')
    #        dtby = pd.Timedelta('1min')
    #        timeby = pd.date_range(btimeby, etimeby, freq=dtby)
    #        epochhour = (timeby - pd.Timestamp('2010-03-23 00:00:00'))\
    #                /pd.Timedelta('1hour')
    #        fnimf1 = '/home/guod/WD4T/gitm/run_imfby/run1c/imf1.dat'
    #        imf1 = pd.read_csv(
    #                fnimf1, delim_whitespace=True, comment='#',
    #                header=None,
    #                names=('year', 'month', 'day', 'hour', 'minute', 'second',
    #                       'ms', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
    #                usecols=['by'])
    #        imf1 = pd.DataFrame(np.array(imf1), index=timeby, columns=['By'])
    #        plt.plot(epochhour, imf1.By, 'b')
    #        fnimf2 = '/home/guod/WD4T/gitm/run_imfby/run2c/imf2.dat'
    #        imf2 = pd.read_csv(
    #                fnimf2, delim_whitespace=True, comment='#',
    #                header=None,
    #                names=('year', 'month', 'day', 'hour', 'minute', 'second',
    #                       'ms', 'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
    #                usecols=['by'])
    #        imf2 = pd.DataFrame(np.array(imf2), index=timeby, columns=['By'])
    #        plt.plot(epochhour, imf2.By, 'r')
    #        plt.xlim(-1, 5)
    #        plt.xticks(np.arange(-1, 5.1, 0.5))
    #        plt.xlabel('Epoch (hour)')
    #        plt.ylim([-10, 10])
    #        plt.ylabel(r'IMF $B_Y$')
    #        plt.axvline(timeidxepoch[i], color='k', linestyle='--')

    #        # read gitm data
    #        g1, g2 = [gitm.GitmBin(
    #                k[i], varlist=['Rho', 'Temperature', 'V!Dn!N (north)',
    #                               'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
    #                               'V!Di!N (north)']) for k in [fn1, fn2]]

    #        # create axis
    #        ax = list(range(6))
    #        projection = ax.copy()
    #        for ins in range(2):
    #            nlat, slat = [90, 50] if ins==0 else [-50, -90]
    #            for irun in range(3):
    #                ax[ins+irun*2], projection[ins+irun*2] = gcc.create_map(
    #                        4, 2, 3+ins+irun*2, 'polar', nlat=nlat, slat=slat,
    #                        dlat=10, centrallon=g3ca.calculate_centrallon(
    #                            g1, 'polar',  useLT=True),
    #                        coastlines=False)

    #        # Temperature
    #        lon1, lat1, zdata1 = g3ca.contour_data('Temperature', g1, alt=400)
    #        lon2, lat2, zdata2 = g3ca.contour_data('Temperature', g2, alt=400)
    #        hc = [ax[k].contourf(
    #                lon1, lat1, zdata1, transform=ccrs.PlateCarree(),
    #                levels=np.linspace(1300, 1550, 21), cmap='viridis',
    #                extend='both') for k in [0, 1]]
    #        hc = [ax[k].contourf(
    #                lon2, lat2, zdata2, transform=ccrs.PlateCarree(),
    #                levels=np.linspace(1300, 1550, 21), cmap='viridis',
    #                extend='both') for k in [2, 3]]
    #        diffzdata = zdata2-zdata1
    #        hc = [ax[k].contourf(
    #                lon2, lat2, diffzdata, transform=ccrs.PlateCarree(),
    #                levels=np.linspace(-200, 200, 21), cmap='seismic',
    #                extend='both') for k in [4, 5]]

    #        # magnetic poles
    #        for iax in range(6):
    #            ax[iax].scatter(
    #                    qlonn, qlatn, color='k', transform=ccrs.PlateCarree())
    #            ax[iax].scatter(
    #                    qlons, qlats, color='k', transform=ccrs.PlateCarree())
    #        return
    #    anim = animation.FuncAnimation(
    #            fig, animate_temperature, interval=200, frames=len(fn1))
    #    anim.save(path+'temperature.mp4', writer='ffmpeg')
    #    plt.close('all')
    return


def snapshot_30min():
    # Every 30 minutes, see the results.
    stime = pd.Timestamp('2010-03-23 00:00:00')
    etime = pd.Timestamp('2010-03-23 01:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fn1 = ['/home/guod/big/raid4/guod/run_imfby/run1c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]
    fn2 = ['/home/guod/big/raid4/guod/run_imfby/run2c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]
    #oo = pd.read_pickle('/home/guod/data/tmp/w3_05_02.dat')
    nlat, slat, alt = 90, 50, 400
    path =  '/home/guod/Documents/work/fig/w03/snapshot/'
    #altarray = oo.index.get_level_values(1).values
    #ialt = np.argmin(np.abs(altarray-alt*1000))
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)
    for k0 in range(len(fn1)):
        # g1 = gitm.GitmBin(
        #         fn1[k0],
        #         varlist=['Rho', 'Temperature',
        #                  'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
        #                  'V!Di!N (east)', 'V!Di!N (north)'])
        g2 = gitm.GitmBin(
                fn2[k0],
                varlist=['Rho', 'Temperature',
                         'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
                         'V!Di!N (east)', 'V!Di!N (north)'])
        gd.calc_divergence(g2)
        # location of density minima
        #oot = oo.loc[(g2['time'], oo.index.get_level_values(1)[ialt]), :]
        #lat, lon = oot['lat'], oot['lon']
        #lt = g2['time'].hour+g2['time'].minute/60+g2['time'].second/3600+lon/15
        #theta, r = lt/12*np.pi, 90-lat

        #    # Density and wind run1
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon0, lat0, zdata0 = g3ca.contour_data('Rho', g1, alt=400)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
        #                     extend='both')
        #    lon0, lat0, ewind, nwind = g3ca.vector_data(g1, 'neutral', alt=400)
        #    lon0, lat0, ewind, nwind = g3ca.convert_vector(
        #            lon0, lat0, ewind, nwind, plot_type='polar',
        #            projection=projection)
        #    hq = ax.quiver(
        #            lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
        #            alpha=0.5, regrid_shape=20)
        #    ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
        #    hc = plt.colorbar(hc, ticks=np.arange(3, 7)*1e-12)
        #    hc.set_label(r'$\rho$ (kg/m$^3$)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_den_win_run1_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # Density and wind run2
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon0, lat0, zdata0 = g3ca.contour_data('Rho', g2, alt=400)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
        #                     extend='both')
        #    lon0, lat0, ewind, nwind = g3ca.vector_data(g2, 'neutral', alt=400)
        #    lon0, lat0, ewind, nwind = g3ca.convert_vector(
        #            lon0, lat0, ewind, nwind, plot_type='polar',
        #            projection=projection)
        #    hq = ax.quiver(
        #            lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
        #            alpha=0.5, regrid_shape=20)
        #    ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
        #    hc = plt.colorbar(hc, ticks=np.arange(3, 7)*1e-12)
        #    hc.set_label(r'$\rho$ (kg/m$^3$)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_den_win_run2_'+
        #            g2['time'].strftime('%H%M')+'_new.pdf')

        #    # Density and wind difference
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=400)
        #    lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2, alt=400)
        #    zdata0 = 100*(zdata2-zdata1)/zdata1
        #    hc = ax.contourf(lon1, lat1, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(-30, 30, 21), cmap='seismic',
        #                     extend='both')
        #    lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=400)
        #    lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=400)
        #    lon0, lat0, ewind, nwind =g3ca.convert_vector(
        #            lon1, lat1, ewind2-ewind1, nwind2-nwind1,
        #            plot_type='polar', projection=projection)
        #    hq = ax.quiver(
        #            lon0, lat0, ewind, nwind, scale=1000, scale_units='inches',
        #            alpha=0.5, regrid_shape=20)
        #    ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
        #    hc = plt.colorbar(hc, ticks = np.arange(-30, 31, 10))
        #    hc.set_label(r'$100\times\frac{\rho_2-\rho_1}{\rho_1}$')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_den_win_diff_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # Temperature run1
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g1, alt=400)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(1300, 1550, 21), cmap='viridis',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(1300, 1551, 50))
        #    hc.set_label('T (degree)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_temperature_run1_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # Temperature run2
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g2, alt=400)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(1300, 1550, 21), cmap='viridis',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(1300, 1551, 50))
        #    hc.set_label('T (degree)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_temperature_run2_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # Temperature run2 - run1
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
        #            coastlines=False)
        #    lon1, lat1, zdata1 = g3ca.contour_data('Temperature', g1, alt=400)
        #    lon2, lat2, zdata2 = g3ca.contour_data('Temperature', g2, alt=400)
        #    zdata0 = zdata2-zdata1
        #    hc = ax.contourf(lon1, lat1, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(-200, 200, 21), cmap='seismic',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(-200, 201, 50))
        #    hc.set_label('T$_2$ - T$_1$ (degree)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_temperature_diff_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # neutral and ion velocity difference, Run 1
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
        #            coastlines=False, useLT=True)
        #    lon0, lat0, nnwind = g3ca.contour_data('V!Dn!N (north)', g1, alt=400)
        #    lon0, lat0, newind = g3ca.contour_data('V!Dn!N (east)', g1, alt=400)
        #    lon0, lat0, inwind = g3ca.contour_data('V!Di!N (north)', g1, alt=400)
        #    lon0, lat0, iewind = g3ca.contour_data('V!Di!N (east)', g1, alt=400)
        #    zdata0 = np.sqrt((nnwind-inwind)**2+(newind-iewind)**2)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(0, 500, 21), cmap='viridis',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(0, 501, 100))
        #    hc.set_label('Vn-i (m/s)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_dv_run1_'+
        #            g1['time'].strftime('%H%M')+'_new.pdf')

        #    # neutral and ion velocity difference, Run 2
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
        #            coastlines=False, useLT=True)
        #    lon0, lat0, nnwind = g3ca.contour_data('V!Dn!N (north)', g2, alt=400)
        #    lon0, lat0, newind = g3ca.contour_data('V!Dn!N (east)', g2, alt=400)
        #    lon0, lat0, inwind = g3ca.contour_data('V!Di!N (north)', g2, alt=400)
        #    lon0, lat0, iewind = g3ca.contour_data('V!Di!N (east)', g2, alt=400)
        #    zdata0 = np.sqrt((nnwind-inwind)**2+(newind-iewind)**2)
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(0, 500, 21), cmap='viridis',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(0, 501, 100))
        #    hc.set_label('Vn-i (m/s)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_dv_run2_'+
        #            g2['time'].strftime('%H%M')+'_new.pdf')

        #    # neutral and ion velocity difference, Run 2 - Run 1
        #    plt.figure()
        #    ax, projection = gcc.create_map(
        #            1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
        #            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
        #            coastlines=False, useLT=True)
        #    lon1, lat1, nnwind1 = g3ca.contour_data('V!Dn!N (north)', g1, alt=400)
        #    lon1, lat1, newind1 = g3ca.contour_data('V!Dn!N (east)', g1, alt=400)
        #    lon1, lat1, inwind1 = g3ca.contour_data('V!Di!N (north)', g1, alt=400)
        #    lon1, lat1, iewind1 = g3ca.contour_data('V!Di!N (east)', g1, alt=400)
        #    zdata1 = np.sqrt((nnwind1-inwind1)**2+(newind1-iewind1)**2)
        #    lon2, lat2, nnwind2 = g3ca.contour_data('V!Dn!N (north)', g2, alt=400)
        #    lon2, lat2, newind2 = g3ca.contour_data('V!Dn!N (east)', g2, alt=400)
        #    lon2, lat2, inwind2 = g3ca.contour_data('V!Di!N (north)', g2, alt=400)
        #    lon2, lat2, iewind2 = g3ca.contour_data('V!Di!N (east)', g2, alt=400)
        #    zdata2 = np.sqrt((nnwind2-inwind2)**2+(newind2-iewind2)**2)
        #    zdata0 = zdata2 -zdata1
        #    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),
        #                     levels=np.linspace(-250, 250, 21), cmap='seismic',
        #                     extend='both')
        #    hc = plt.colorbar(hc, ticks=np.arange(-250, 251, 50))
        #    hc.set_label(r'$\Delta$Vn-i (m/s)')
        #    ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        #    plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        #    plt.savefig(
        #            path+'w03_func08_dv_diff_'+
        #            g2['time'].strftime('%H%M')+'_new.pdf')

        # advection term
        ggt.calc_gradient(g2, 'Rho', 'gradn', component='north')
        ggt.calc_gradient(g2, 'Rho', 'grade', component='east')
        ggt.calc_gradient(g2, 'Rho', 'gradr', component='radial')
        gd.calc_divergence(g2, neuion='neutral', name='divergence')

        lon0, lat0, ewind0 = g3ca.contour_data('V!Dn!N (east)', g2, alt=alt)
        lon0, lat0, nwind0 = g3ca.contour_data('V!Dn!N (north)', g2, alt=alt)
        lon0, lat0, uwind0 = g3ca.contour_data('V!Dn!N (up)', g2, alt=alt)
        ewind0 = (ewind0 +
                  ((2*np.pi)/(24*3600))*(Re+alt*1000)+np.cos(lat0*np.pi/180))
        lon0, lat0, gradn0 = g3ca.contour_data('gradn', g2, alt=alt)
        lon0, lat0, grade0 = g3ca.contour_data('grade', g2, alt=alt)
        lon0, lat0, gradr0 = g3ca.contour_data('gradr', g2, alt=alt)

        lon0, lat0, rho0 = g3ca.contour_data('Rho', g2, alt=alt)
        lon0, lat0, div0 = g3ca.contour_data('divergence', g2, alt=alt)

        hadvect = ewind0 * grade0 + nwind0 * gradn0
        vadvect = uwind0*gradr0

        eee = rho0*div0

        plt.figure()
        ax, projection = gcc.create_map(
                1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -(hadvect+vadvect+eee),
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-8e-15, 8e-15, 21), cmap='seismic')
        hcb = plt.colorbar(hc, ticks=np.arange(-8, 9, 2)*1e-15)
        hcb.set_label(r'$-\nabla\cdot(\rho \vec{v})$')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        plt.savefig(path+'w03_func08_den_change_all_run2_'+
                    g2['time'].strftime('%H%M')+'_new.pdf')

        plt.figure()
        ax, projection = gcc.create_map(
                1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -(hadvect+vadvect),
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-8e-15, 8e-15, 21), cmap='seismic')
        hcb = plt.colorbar(hc, ticks=np.arange(-8, 9, 2)*1e-15)
        hcb.set_label(r'$-\vec{v}\cdot\nabla\rho$')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        plt.savefig(path+'w03_func08_den_change_advect_run2_'+
                    g2['time'].strftime('%H%M')+'_new.pdf')

        plt.figure()
        ax, projection = gcc.create_map(
                1, 1, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g2, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -eee,
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-8e-15, 8e-15, 21), cmap='seismic')
        hcb = plt.colorbar(hc, ticks=np.arange(-8, 9, 2)*1e-15)
        hcb.set_label(r'$-\rho\nabla\cdot\vec{v}$')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title(g2['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        plt.savefig(path+'w03_func08_den_change_volume_run2_'+
                    g2['time'].strftime('%H%M')+'_new.pdf')

        plt.close('all')
    return


def initial_condition():
    fn = ('/home/guod/big/raid4/guod/run_imfby/run1c/data/'
          '3DALL_t100323_000000.bin')
    savepath =  '/home/guod/Documents/work/work3paper/figs/'
    nlat, slat, alt = 90, 50, 400

    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    g = gitm.GitmBin(
            fn, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])

    fig = plt.figure(figsize=(8, 5))
    # density and wind
    ax, projection = gcc.create_map(
            1, 2, 1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
            coastlines=False)
    lon0, lat0, rho0 = g3ca.contour_data('Rho', g, alt=alt)
    hc = ax.contourf(lon0, lat0, rho0, transform=ccrs.PlateCarree(),
                     levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
                     extend='both')
    lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
    lon0, lat0, ewind, nwind = g3ca.convert_vector(
            lon0, lat0, ewind, nwind, plot_type='polar',
            projection=projection)
    hq = ax.quiver(
            lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
            alpha=0.5, regrid_shape=20)
    ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
    cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=-0.01)
    hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(3, 7)*1e-12,
                       orientation='horizontal')
    xticklb = [r'$%d\times10^{-12}$' % k for k in np.arange(3, 7)]
    hcb.ax.set_xticklabels(xticklb)
    hcb.set_label(r'$\rho$ (kg/m$^3$)')
    # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
    ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    ax.text(0, 1.05, '( a )', transform=ax.transAxes)

    # temperature
    ax, projection = gcc.create_map(
            1, 2, 2, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
            coastlines=False)
    lon0, lat0, temp0 = g3ca.contour_data('Temperature', g, alt=alt)
    hc = ax.contourf(lon0, lat0, temp0, transform=ccrs.PlateCarree(),
                     levels=np.linspace(1300, 1550, 21), cmap='viridis',
                     extend='both')
    cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=-0.01)
    hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(1300, 1551, 50),
                       orientation='horizontal')
    hcb.set_label('T (K)')
    # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
    ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    ax.text(0, 1.05, '( b )', transform=ax.transAxes)

    plt.savefig(savepath+'figure1.eps')

    return


def half_hour_condition():
    fn1 = ('/home/guod/big/raid4/guod/run_imfby/run1c/data/'
          '3DALL_t100323_003000.bin')
    fn2 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_003000.bin')
    savepath =  '/home/guod/Documents/work/work3paper/figs/'
    nlat, slat, alt = 90, 50, 400

    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    g1 = gitm.GitmBin(
            fn1, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])
    g2 = gitm.GitmBin(
            fn2, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])

    fig = plt.figure(figsize=(9.3, 8.2))

    # run 1 and run 2
    for k11, g in enumerate([g1, g2]):
        if k11==0:
            abcd = '( a )'
        else:
            abcd = '( b )'
        ax, projection = gcc.create_map(
                2, 3, k11+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, rho0 = g3ca.contour_data('Rho', g, alt=alt)
        hc = ax.contourf(lon0, lat0, rho0, transform=ccrs.PlateCarree(),
                         levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
                         extend='both')
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                alpha=0.5, regrid_shape=20)
        ax.quiverkey(hq, 0.93, -0.05, 1000, '1000 m/s')
        cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.0)
        hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(3, 7)*1e-12,
                           orientation='horizontal')
        xticklb = [r'$%d\times10^{-12}$' % k for k in np.arange(3, 7)]
        hcb.ax.set_xticklabels(xticklb)
        if k11==0:
            hcb.set_label(r'$\rho_1$ (kg/m$^3$)')
        else:
            hcb.set_label(r'$\rho_2$ (kg/m$^3$)')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1.05, abcd, transform=ax.transAxes)

    for k11, g in enumerate([g1, g2]):
        if k11==0:
            abcd = '( d )'
        else:
            abcd = '( e )'
        ax, projection = gcc.create_map(
                2, 3, k11+4, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, temp0 = g3ca.contour_data('Temperature', g, alt=alt)
        hc = ax.contourf(lon0, lat0, temp0, transform=ccrs.PlateCarree(),
                         levels=np.linspace(1300, 1550, 21), cmap='viridis',
                         extend='both')
        cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.0)
        hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(1300, 1551, 50),
                           orientation='horizontal')
        if k11==0:
            hcb.set_label(r'$T_1$ (K)')
        else:
            hcb.set_label(r'$T_2$ (K)')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1.05, abcd, transform=ax.transAxes)

    ax, projection = gcc.create_map(
            2, 3, 3, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
            coastlines=False)
    lon1, lat1, rho1 = g3ca.contour_data('Rho', g1, alt=alt)
    lon2, lat2, rho2 = g3ca.contour_data('Rho', g2, alt=alt)
    lon0, lat0, rho0 = lon1, lat1, 100*(rho2-rho1)/rho1
    hc = ax.contourf(lon0, lat0, rho0, transform=ccrs.PlateCarree(),
                     levels=np.linspace(-30, 30, 21), cmap='seismic',
                     extend='both')
    lon0, lat0, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=alt)
    lon0, lat0, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=alt)
    ewind0, nwind0 = ewind2-ewind1, nwind2-nwind1
    lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind0, nwind0, plot_type='polar',
            projection=projection)
    hq = ax.quiver(
            lon0, lat0, ewind0, nwind0, scale=1500, scale_units='inches',
            alpha=0.5, regrid_shape=20)
    ax.quiverkey(hq, 0.93, -0.05, 1000, '1000 m/s')
    cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.0)
    hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(-30, 30, 10),
                       orientation='horizontal')
    hcb.set_label(r'$100\times\frac{\rho_2-\rho_1}{\rho_1}$')
    # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
    ax.set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    ax.text(0, 1.05, '( c )', transform=ax.transAxes)

    ax, projection = gcc.create_map(
            2, 3, 6, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
            coastlines=False)
    lon1, lat1, temp1 = g3ca.contour_data('Temperature', g1, alt=alt)
    lon2, lat2, temp2 = g3ca.contour_data('Temperature', g2, alt=alt)
    lon0, lat0, temp0 = lon1, lat1, temp2-temp1
    hc = ax.contourf(lon0, lat0, temp0, transform=ccrs.PlateCarree(),
                     levels=np.linspace(-200, 200, 21), cmap='seismic',
                     extend='both')
    cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.0)
    hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(-200, 201, 100),
                       orientation='horizontal')
    hcb.set_label(r'$T_2-T_1$ (K)')
    # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
    ax.set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    ax.text(0, 1.05, '( f )', transform=ax.transAxes)

    plt.savefig(savepath+'figure2.eps')
    return


def time_evolution_mean_density_temperature():
    # UT change in density and temperature
    dt = pd.date_range(
            start='2010-03-23 00:00:00', end='2010-03-23 05:00:00',
            freq='5min')
    nlat, slat, alt = 90, 50, 400
    if False: # data preparation
        whichrun = 'run2'
        path = '/home/guod/big/raid4/guod/run_imfby/{:s}c/data/'.format(
                whichrun)
        fn = (path+'3DALL_t{:s}.bin'.format(k.strftime('%y%m%d_%H%M%S'))
              for k in dt)
        oo = pd.DataFrame(columns=('time', 'rhomean', 'tempmean', 'vdiffmean',
                                   'itempmean'))
        for k00, k0 in enumerate(fn):
            if not os.path.isfile(k0):
                continue
            g = gitm.GitmBin(
                    k0,
                    varlist=['Rho', 'Temperature', 'V!Dn!N (north)',
                             'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                             'V!Di!N (north)', 'iTemperature'])
            ilat = (g['dLat'][0, :, 0]>slat) & (g['dLat'][0, :, 0]<nlat)
            ilon = (g['dLon'][:, 0, 0]>=0) & (g['dLon'][:, 0, 0]<=360)
            ialt = np.argmin(np.abs(g['Altitude'][0, 0, :]-alt*1000))
            gt = g['time']
            rho, temp, lat, itemp = (
                    g[k][ilon][:, ilat][:, :, ialt]
                    for k in ('Rho', 'Temperature',
                              'Latitude', 'iTemperature'))
            nnwind, newind, inwind, iewind = (
                    g[k][ilon][:, ilat][:, :, ialt]
                    for k in ( 'V!Dn!N (up)', 'V!Dn!N (east)',
                               'V!Di!N (north)', 'V!Di!N (east)'))
            vdiff = np.sqrt((nnwind-inwind)**2 + (newind-iewind)**2)
            vdiffmean = np.mean(vdiff*np.cos(lat))/np.mean(np.cos(lat))
            rhomean = np.mean(rho*np.cos(lat))/np.mean(np.cos(lat))
            tempmean = np.mean(temp*np.cos(lat))/np.mean(np.cos(lat))
            itempmean = np.mean(itemp*np.cos(lat))/np.mean(np.cos(lat))
            oo = oo.append(pd.DataFrame(
                [[gt, rhomean, tempmean, vdiffmean, itempmean]],
                columns=('time', 'rhomean', 'tempmean',
                         'vdiffmean', 'itempmean')))
        pd.to_pickle(
                oo, '/home/guod/data/tmp/w3_09_{:s}.dat'.format(whichrun[-1]))
    oo1 = pd.read_pickle('/home/guod/data/tmp/w3_09_1.dat')
    oo2 = pd.read_pickle('/home/guod/data/tmp/w3_09_2.dat')
    xhour = pd.TimedeltaIndex(
            oo1['time']-pd.Timestamp('2010-03-23 00:00:00')).\
            total_seconds()/3600

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(7, 5.5))
    plt.sca(ax[0])
    plt.plot(xhour, oo1['rhomean'], label='Run 1')
    plt.plot(xhour, oo2['rhomean'], label='Run 2')
    plt.legend()
    plt.xlim(0, 5)
    plt.ylim(5.05e-12, 5.2e-12)

    plt.sca(ax[1])
    plt.plot(xhour, oo1['tempmean'])
    plt.plot(xhour, oo2['tempmean'])
    plt.xlim(0, 5)
    plt.ylim(1320, 1380)
    plt.xticks(np.arange(0, 4.1, 0.5))

    plt.sca(ax[2])
    plt.plot(xhour, oo1['itempmean'])
    plt.plot(xhour, oo2['itempmean'])
    plt.xlim(0, 5)
    plt.ylim(1400, 1525)
    plt.xticks(np.arange(0, 5.1, 0.5))
    mf.subplots_xylabel(
            ax, 3, 1, xlabel='Epoch (hour)',
            ylabel=[r'$\rho$ $(kg/m^3)$', '$T_n$ (K)', '$T_i$ (K)'],
            fontweight='bold')
    mf.subplots_create_abcde(ax, direction='row', x=0.95, y=1.05)
    plt.savefig('/home/guod/Documents/work/work3paper/figs/figure3.eps')
    return


def density_change_reason():
    savepath =  '/home/guod/Documents/work/work3paper/figs/'
    nlat, slat, alt = 90, 50, 400

    fn1 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_000500.bin')
    fn2 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_002500.bin')

    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=alt)

    g1 = gitm.GitmBin(
            fn1, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])
    g2 = gitm.GitmBin(
            fn2, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])

    fig = plt.figure(figsize=(9.3, 6.36))
    abcd = ['( %s )' % k for k in 'abcdef']
    for k11, g in enumerate([g1, g2]):
        ggt.calc_gradient(g, 'Rho', 'gradn', component='north')
        ggt.calc_gradient(g, 'Rho', 'grade', component='east')
        ggt.calc_gradient(g, 'Rho', 'gradr', component='radial')
        gd.calc_divergence(g, neuion='neutral', name='divergence')

        lon0, lat0, ewind0 = g3ca.contour_data('V!Dn!N (east)', g, alt=alt)
        lon0, lat0, nwind0 = g3ca.contour_data('V!Dn!N (north)', g, alt=alt)
        lon0, lat0, uwind0 = g3ca.contour_data('V!Dn!N (up)', g, alt=alt)
        ewind0 = (ewind0 +
                  ((2*np.pi)/(24*3600))*(Re+alt*1000)+np.cos(lat0*np.pi/180))
        lon0, lat0, gradn0 = g3ca.contour_data('gradn', g, alt=alt)
        lon0, lat0, grade0 = g3ca.contour_data('grade', g, alt=alt)
        lon0, lat0, gradr0 = g3ca.contour_data('gradr', g, alt=alt)

        lon0, lat0, rho0 = g3ca.contour_data('Rho', g, alt=alt)
        lon0, lat0, div0 = g3ca.contour_data('divergence', g, alt=alt)

        advect = ewind0 * grade0 + nwind0 * gradn0 + uwind0 * gradr0

        volume_change = rho0*div0

        ax, projection = gcc.create_map(
                2, 3, 1+k11*3, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -advect,
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-2e-15, 2e-15, 21), cmap='seismic')
        if k11==1:
            cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.09)
            hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(-8, 9, 1)*1e-15,
                               orientation='horizontal')
            hcb.set_label(r'$\vec{v}\cdot\nabla\rho$')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1, abcd[0+k11*3], transform=ax.transAxes)

        ax, projection = gcc.create_map(
                2, 3, 2+k11*3, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -volume_change,
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-2e-15, 2e-15, 21), cmap='seismic')
        if k11==1:
            cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.09)
            hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(-8, 9, 1)*1e-15,
                               orientation='horizontal')
            hcb.set_label(r'$\rho\nabla\cdot\vec{v}$')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1, abcd[1+k11*3], transform=ax.transAxes)

        ax, projection = gcc.create_map(
                2, 3, 3+k11*3, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        hc = ax.contourf(
                lon0, lat0, -(advect+volume_change),
                transform=ccrs.PlateCarree(), extend='both',
                levels = np.linspace(-2e-15, 2e-15, 21), cmap='seismic')
        if k11==1:
            cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=0.09)
            hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(-8, 9, 1)*1e-15,
                               orientation='horizontal')
            hcb.set_label(r'$-\nabla\cdot(\rho\vec{v})$')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1, abcd[2+k11*3], transform=ax.transAxes)

        plt.savefig(savepath+'figure4.eps')
    return


def min_density_location():
    fn1 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_000000.bin')
    fn2 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_013500.bin')
    fn3 = ('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
          '3DALL_t100323_030000.bin')
    savepath =  '/home/guod/Documents/work/work3paper/figs/'
    nlat, slat, alt = 90, 50, 400

    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    g1 = gitm.GitmBin(
            fn1, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])
    g2 = gitm.GitmBin(
            fn2, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])
    g3 = gitm.GitmBin(
            fn3, varlist=[
                'Rho', 'Temperature', 'V!Dn!N (north)',
                'V!Dn!N (east)', 'V!Dn!N (up)', 'V!Di!N (east)',
                'V!Di!N (north)'])

    fig = plt.figure(figsize=(9.3, 4.2))

    # run 1 and run 2
    for k11, g in enumerate([g1, g2, g3]):
        if k11==0:
            abcd = '( a )'
        elif k11==1:
            abcd = '( b )'
        else:
            abcd = '( c )'
        ax, projection = gcc.create_map(
                1, 3, k11+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, rho0 = g3ca.contour_data('Rho', g, alt=alt)
        hc = ax.contourf(lon0, lat0, rho0, transform=ccrs.PlateCarree(),
                         levels=np.linspace(3e-12, 6e-12, 21), cmap='viridis',
                         extend='both')
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                alpha=0.5, regrid_shape=20)
        ax.quiverkey(hq, 0.93, -0.05, 1000, '1000 m/s')
        cax = mf.subplots_create_cbar_axis(ax, 'bottom', pad=-0.01)
        hcb = plt.colorbar(hc, cax=cax, ticks=np.arange(3, 7)*1e-12,
                           orientation='horizontal')
        xticklb = [r'$%d\times10^{-12}$' % k for k in np.arange(3, 7)]
        hcb.ax.set_xticklabels(xticklb)
        hcb.set_label(r'$\rho_2$ (kg/m$^3$)')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.set_title(g['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax.text(0, 1.05, abcd, transform=ax.transAxes)
    plt.savefig(savepath+'figure5.eps')
    return


def time_evolution_min_density():
    sns.set_style({'ytick.major.size':3, 'xtick.major.size':3})
    dt = pd.date_range(
            start='2010-03-23 00:00:00', end='2010-03-23 05:00:00',
            freq='10min')
    dthour = (dt - pd.Timestamp('2010-03-23 00:00:00'))/pd.Timedelta('1hour')
    nlat, slat = 90, 50
    if False: # data preparation
        whichrun = 'run1'
        path = '/home/guod/big/raid4/guod/run_imfby/{:s}c/data/'.format(
                whichrun)
        fn = (path+'3DALL_t{:s}.bin'.format(k.strftime('%y%m%d_%H%M%S'))
              for k in dt)
        oorho = []
        oolat = []
        for k00, k0 in enumerate(fn):
            if not os.path.isfile(k0):
                continue
            g = gitm.GitmBin(
                    k0,
                    varlist=['Rho', 'Temperature', 'iTemperature',
                             'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
                             'V!Di!N (east)', 'V!Di!N (north)'])
            ilat = (g['dLat'][0, :, 0]>slat) & (g['dLat'][0, :, 0]<nlat)
            ilon = (g['dLon'][:, 0, 0]>=0) & (g['dLon'][:, 0, 0]<=360)
            gt = g['time']
            rho = g['Rho'][ilon, ...][:, ilat, :]
            lat = g['Latitude'][ilon, ...][:, ilat, :]
            lon = g['Longitude'][ilon, ...][:, ilat, :]
            alt = g['Altitude'][0, 0, :]
            oorhot = alt.copy()
            oolatt = alt.copy()
            for k11, k1 in enumerate(alt):
                rhominidx = np.unravel_index(
                        np.argmin(rho[..., k11]), rho[..., k11].shape)
                rhomin = np.min(rho[..., k11])
                latmin = (lat[rhominidx[0], rhominidx[1], k11])/np.pi*180
                lonleftidx = (rhominidx[0]-11)%lon.shape[0]
                lonrightidx = (rhominidx[0]+11)%lon.shape[0]
                rhomean = (rho[lonleftidx, rhominidx[1], k11] +
                           rho[lonrightidx, rhominidx[1], k11])/2
                oorhot[k11] = rhomin/rhomean
                oolatt[k11] = latmin
            oorhot = oorhot.reshape(-1, 1)
            oorho.append(oorhot)
            oolatt = oolatt.reshape(-1, 1)
            oolat.append(oolatt)
        oorho = np.concatenate(oorho, axis=1)
        oolat = np.concatenate(oolat, axis=1)
        pd.to_pickle(
                [dthour, alt, oorho, oolat],
                '/home/guod/data/tmp/w3_min_density_{:s}.dat'.\
                format(whichrun[-1]))
    dthour1, alt1, oorho1, oolat1 = pd.read_pickle(
            '/home/guod/data/tmp/w3_min_density_1.dat')
    dthour2, alt2, oorho2, oolat2 = pd.read_pickle(
            '/home/guod/data/tmp/w3_min_density_2.dat')

    fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8.3, 6))

    plt.sca(ax[0, 0])
    hc1 = plt.contourf(dthour1, alt1/1000, oorho1, cmap='viridis',
                 levels=np.linspace(0.5, 1,  21))
    plt.title('Run 1', fontweight='bold', fontsize='large')

    plt.sca(ax[0, 1])
    hc2 = plt.contourf(dthour1, alt1/1000, oorho2, cmap='viridis',
                 levels=np.linspace(0.5, 1,  21))
    plt.title('Run 2', fontweight='bold', fontsize='large')

    plt.sca(ax[1, 0])
    hc3 = plt.contourf(dthour1, alt1/1000, oolat1, cmap='viridis',
                 levels=np.linspace(50, 85,  21))

    plt.sca(ax[1, 1])
    hc4 = plt.contourf(dthour1, alt1/1000, oolat2, cmap='viridis',
                 levels=np.linspace(50, 85,  21))
    plt.xlim(0, 5)
    plt.ylim(100, 700)

    mf.subplots_xylabel(ax, 2, 2, xlabel='Epoch (hour)',
                        ylabel='Altitude (km)', fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.90, 1])
    cax1 = mf.subplots_create_cbar_axis(ax[0, 1])
    plt.colorbar(hc2, cax=cax1, ticks=np.arange(0.5, 1.1, 0.1))
    plt.ylabel(r'$\rho_{min}/\rho_b$', fontsize='large', fontweight='bold')
    cax2 = mf.subplots_create_cbar_axis(ax[1, 1])
    plt.colorbar(hc4, cax=cax2, ticks=np.arange(50, 90, 5))
    plt.ylabel('Latitude', fontsize='large', fontweight='bold')

    mf.subplots_create_abcde(ax, direction='row', y=1.02)
    plt.savefig('/home/guod/Documents/work/work3paper/figs/figure6.eps')
    return
# END
#------------------------------------------------------------
if __name__ == '__main__':
    plt.close('all')
    a = density_change_reason()
    plt.show()
