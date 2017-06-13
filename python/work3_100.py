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
from apexpy import Apex
from matplotlib.ticker import AutoMinorLocator
sns.set('paper', 'whitegrid')

DATADIR = os.environ.get('DATAPATH')
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
        fig1.savefig(figdir + 'w03_func03_03a')

    lb = ['CHAMP', 'GRACE', 'GOCE']
    cl = list('rbk')
    if True:  # Read satellite data
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


def func5():
    # Analyze gitm output from run_imfby
    import os
    import gitm
    import gitm_3D_const_alt as g3ca
    import matplotlib.animation as animation
    import gitm_vorticity as gv
    from apexpy import Apex
    from matplotlib.ticker import AutoMinorLocator

    # Find file names (100322_230000 -> 100323_060000)
    bdir = '/home/guod/WD4T/gitm/run_imfby/run1b/data/'
    tname = os.listdir(bdir)
    tname = [i for i in tname if '.bin' in os.path.splitext(i)]
    tname = np.sort(tname)
    inds, = np.where(tname=='3DALL_t100322_230000.bin')[0]
    #inds, = np.where(tname=='3DALL_t100322_000000.bin')[0]
    basename = tname[inds:]

    pdir = '/home/guod/WD4T/gitm/run_imfby/run2b/data/'
    tname = os.listdir(pdir)
    tname = [i for i in tname if '.bin' in os.path.splitext(i)]
    tname = np.sort(tname)
    inds, = np.where(tname=='3DALL_t100322_230000.bin')[0]
    pertname = tname[inds:]

    if False:  # Density and wind change in run1 and run2
        zlim = [[1e-9, 2e-9], [1.5e-10, 3e-10], [1.8e-11, 4e-11],
                [0.4e-11, 1e-11], [1.1e-12, 2.7e-12]]
        for kk in ['Run1', 'Run2']:
        #for kk in ['Run1']:#, 'Run2']:
            if kk is 'Run1':
                bpname, bpdir = basename, bdir
            else:
                bpname, bpdir = pertname, pdir
            for k0 in bpname:
                gg = gitm.GitmBin(
                        bpdir+k0,
                        varlist=['Rho', 'V!Dn!N (north)', 'V!Dn!N (east)'])
                for k22, k2 in enumerate([150, 200, 300, 400, 500]):
                    alt = gg['Altitude'][0, 0, :]
                    ialt = np.argmin(abs(alt-k2*1000)) # in GITM, unit of alt is m
                    altx = alt[ialt]/1000
                    # geomagnetic poles
                    mm = Apex(date=2010.3)
                    mnplat, mnplon = mm.convert(
                            90, 0, source='qd', dest='geo', height=altx)
                    mnplt = (
                            gg['time'].hour+gg['time'].minute/60+
                            gg['time'].second/3600)+mnplon/15
                    msplat, msplon = mm.convert(-90, 0, source='qd', dest='geo',
                                                height=altx)
                    msplt = (
                            gg['time'].hour+gg['time'].minute/60+
                            gg['time'].second/3600)+msplon/15
                    for k33, k3 in enumerate(['North', 'South']):
                        fig = plt.figure()
                        ax = plt.subplot(polar=True)
                        print('Figure info: {:s}, {:s}'
                              ' at {:5.1f} km'.format(kk, k3, altx))
                        nlat = 90 if k3 is 'North' else -40
                        slat = -90 if k3 is 'South' else 40
                        rp = 90-mnplat if k3 is 'North' else 90+msplat
                        thetap = (mnplt*np.pi/12 if k3 is 'North'
                                  else (24-msplt)*np.pi/12)
                        ax, hc = g3ca.contour_single(
                                ax, 'Rho', 'pol', gg, alt=k2,
                                nlat=nlat, slat=slat, dlat=10,
                                dlonlt=6, zmax=zlim[k22][1], nzlevels=21,
                                zmin=zlim[k22][0], data_type="contour")
                        ax, hq = g3ca.vector_single(
                                ax, gg, 'neu', 'pol', alt=k2, nlat=nlat,
                                slat=slat, dlat=10, dlonlt=6,
                                color='k', alpha=0.6, scale=1000,
                                scale_units='inches', headwidth=5)
                        ax.quiverkey(hq, 0, 0, 1000, '1000 m/s')
                        ax.scatter(thetap, rp, s=20, c='k')
                        # Time
                        tt = gg['time']
                        ht = plt.title(tt.strftime('%y-%m-%d %H:%M'))
                        ht.set_position(
                                [ht.get_position()[0],
                                 ht.get_position()[1]+0.01])
                        # north or south
                        plt.annotate('{:s}, {:5.1f} km'.format(k3, altx),
                                     [0, 1.09], xycoords='axes fraction',
                                     horizontalalignment='center',
                                     verticalalignment='center')
                        # colorbar
                        hcb = plt.colorbar(hc)#, ticks=np.arange(-50, 51, 10))
                        plt.tight_layout()
                        fig.savefig('/home/guod/WD4T/work3/'+
                                '{:s}_{:s}_{:d}_'.format(kk, k3, k2)+k0+'_den_wind.pdf')
                        plt.close(fig)


    if False:  # Density and wind difference between run2 and run1
        for k0, k1 in zip(basename, pertname):
            g1 = gitm.GitmBin(
                    bdir+k0,
                    varlist=['Rho', 'V!Dn!N (north)', 'V!Dn!N (east)'])
            g2 = gitm.GitmBin(
                    pdir+k1,
                    varlist=['Rho', 'V!Dn!N (north)', 'V!Dn!N (east)'])
            for k22, k2 in enumerate([150, 200, 300, 400, 500]):
                alt = g1['Altitude'][0, 0, :]
                ialt = np.argmin(abs(alt-k2*1000)) # in GITM, unit of alt is m
                altx = alt[ialt]/1000
                # geomagnetic poles
                mm = Apex(date=2010.3)
                mnplat, mnplon = mm.convert(
                        90, 0, source='qd', dest='geo', height=altx)
                mnplt = (
                        g1['time'].hour+
                        g1['time'].minute/60+
                        g1['time'].second/3600)+mnplon/15
                msplat, msplon = mm.convert(
                        -90, 0, source='qd', dest='geo', height=altx)
                msplt = (
                        g1['time'].hour+g1['time'].minute/60+
                        g1['time'].second/3600)+msplon/15
                for k33, k3 in enumerate(['North', 'South']):
                    fig = plt.figure()
                    ax = plt.subplot(polar=True)
                    print('Figure info: {:s} at {:5.1f} km'.format(k3, altx))
                    nlat = 90 if k3 is 'North' else -40
                    slat = -90 if k3 is 'South' else 40
                    rp = 90-mnplat if k3 is 'North' else 90+msplat
                    thetap = (mnplt*np.pi/12
                              if k3 is 'North' else (24-msplt)*np.pi/12)
                    ax, hc = g3ca.contour_diff(
                            ax, 'Rho', 'pol', g1, g2, alt=k2, nzlevels=20,
                            nlat=nlat, slat=slat, dlonlt=6, zmax=40, zmin=-40)
                    ax, hq = g3ca.vector_diff(
                            ax, g1, g2, 'neu', 'pol', alt=k2, nlat=nlat,
                            slat=slat, dlonlt=6, color='k', alpha=0.6,
                            scale=1000, scale_units='inches', headwidth=5)
                    ax.quiverkey(hq, 0, 0, 1000, '1000 m/s')
                    ax.scatter(thetap, rp, s=20, c='k')
                    # Time
                    tt = g2['time']
                    ht = plt.title(tt.strftime('%y-%m-%d %H:%M'))
                    ht.set_position(
                            [ht.get_position()[0], ht.get_position()[1]+0.01])
                    # north or south
                    plt.annotate(
                            '{:s}, {:5.1f} km'.format(k3, altx),
                            [0, 1.09], xycoords='axes fraction',
                            horizontalalignment='center',
                            verticalalignment='center')
                    # colorbar
                    hcb = plt.colorbar(hc, ticks=np.arange(-40, 41, 10))
                    plt.tight_layout()
                    fig.savefig(
                            '/home/guod/WD4T/work3/'+
                            '{:s}_{:d}_'.format(k3, k2)+k1+'_den_wind_diff.pdf')
                    plt.close(fig)

    if False:  # Temperature change in run1 and run2
        zlim = [[650, 850],[1000, 1300], [1100, 1600],
                [1100, 1700], [1100, 1700]]
        for kk in ['Run1', 'Run2']:
        #for kk in ['Run1']:#, 'Run2']:
            if kk is 'Run1':
                bpname, bpdir = basename, bdir
            else:
                bpname, bpdir = pertname, pdir
            for k0 in bpname:
                gg = gitm.GitmBin(bpdir+k0, varlist=['Temperature', 'Rho'])
                for k22, k2 in enumerate([150, 200, 300, 400, 500]):
                    alt = gg['Altitude'][0, 0, :]
                    ialt = np.argmin(abs(alt-k2*1000)) # in GITM, unit of alt is m
                    altx = alt[ialt]/1000
                    # geomagnetic poles
                    mm = Apex(date=2010.3)
                    mnplat, mnplon = mm.convert(
                            90, 0, source='qd', dest='geo', height=altx)
                    mnplt = (gg['time'].hour+ gg['time'].minute/60+
                             gg['time'].second/3600)+mnplon/15
                    msplat, msplon = mm.convert(
                            -90, 0, source='qd', dest='geo', height=altx)
                    msplt = (gg['time'].hour+ gg['time'].minute/60+
                             gg['time'].second/3600)+msplon/15
                    for k33, k3 in enumerate(['North', 'South']):
                        fig = plt.figure()
                        ax = plt.subplot(polar=True)
                        print('Figure info: {:s}, {:s}'
                              ' at {:5.1f} km'.format(kk, k3, altx))
                        nlat = 90 if k3 is 'North' else -40
                        slat = -90 if k3 is 'South' else 40
                        rp = 90-mnplat if k3 is 'North' else 90+msplat
                        thetap = (mnplt*np.pi/12 if k3 is 'North'
                                  else (24-msplt)*np.pi/12)
                        ax, hc = g3ca.contour_single(
                                ax, 'Temperature', 'pol', gg, alt=k2,
                                nlat=nlat, slat=slat, dlonlt=6, nzlevels=20,
                                zmax=zlim[k22][1], zmin=zlim[k22][0])
                        ax.scatter(thetap, rp, s=20, c='k')
                        # Time
                        tt = gg['time']
                        ht = plt.title(tt.strftime('%y-%m-%d %H:%M'))
                        ht.set_position(
                                [ht.get_position()[0],
                                 ht.get_position()[1]+0.01])
                        # north or south
                        plt.annotate('{:s}, {:5.1f} km'.format(k3, altx),
                                     [0, 1.09], xycoords='axes fraction',
                                     horizontalalignment='center',
                                     verticalalignment='center')
                        # colorbar
                        hcb = plt.colorbar(hc)#, ticks=np.arange(-50, 51, 10))
                        plt.tight_layout()
                        fig.savefig('/home/guod/WD4T/work3/'+
                                '{:s}_{:s}_{:d}_'.format(kk, k3, k2)+k0+
                                '_Temperature.pdf')
                        plt.close(fig)

    if False:  # Rho(filled contour) and vorticity change in run1 and run2
        zlim = [[1e-9, 2e-9], [1.5e-10, 3e-10], [1.8e-11, 4e-11],
                [0.4e-11, 1e-11], [1.1e-12, 2.7e-12]]
        for kk in ['Run1', 'Run2']:
        #for kk in ['Run1']:#, 'Run2']:
            if kk is 'Run1':
                bpname, bpdir = basename, bdir
            else:
                bpname, bpdir = pertname, pdir
            for k0 in bpname:
                gg = gitm.GitmBin(
                        bpdir+k0,
                        varlist=['Rho', 'V!Dn!N (north)', 'V!Dn!N (east)'])
                gv.calc_vorticity(gg, 'neutral')
                for k22, k2 in enumerate([150, 200, 300, 400, 500]):
                    alt = gg['Altitude'][0, 0, :]
                    ialt = np.argmin(abs(alt-k2*1000)) # in GITM, unit of alt is m
                    altx = alt[ialt]/1000
                    # geomagnetic poles
                    mm = Apex(date=2010.3)
                    mnplat, mnplon = mm.convert(
                            90, 0, source='qd', dest='geo', height=altx)
                    mnplt = (gg['time'].hour+ gg['time'].minute/60+
                             gg['time'].second/3600)+mnplon/15
                    msplat, msplon = mm.convert(
                            -90, 0, source='qd', dest='geo', height=altx)
                    msplt = (gg['time'].hour+ gg['time'].minute/60+
                             gg['time'].second/3600)+msplon/15
                    for k33, k3 in enumerate(['North', 'South']):
                        fig = plt.figure()
                        ax = plt.subplot(polar=True)
                        print('Figure info: {:s}, {:s}'
                              ' at {:5.1f} km'.format(kk, k3, altx))
                        nlat = 90 if k3 is 'North' else -40
                        slat = -90 if k3 is 'South' else 40
                        rp = 90-mnplat if k3 is 'North' else 90+msplat
                        thetap = (mnplt*np.pi/12 if k3 is 'North'
                                  else (24-msplt)*np.pi/12)
                        ax, hc = g3ca.contour_single(
                                ax, 'Rho', 'pol', gg, alt=k2,
                                nlat=nlat, slat=slat, dlonlt=6, nzlevels=20,
                                zmax=zlim[k22][1], zmin=zlim[k22][0])
                        ax, hv = g3ca.contour_single(
                                ax, 'nvorticity', 'pol', gg, alt=k2, nlat=nlat,
                                slat=slat, dlonlt=6, linewidths=1, alpha=0.8,
                                nzlevels=10, colors='k', fill=False)
                        ax.scatter(thetap, rp, s=20, c='k')
                        # Time
                        tt = gg['time']
                        ht = plt.title(tt.strftime('%y-%m-%d %H:%M'))
                        ht.set_position(
                                [ht.get_position()[0],
                                 ht.get_position()[1]+0.01])
                        # north or south
                        plt.annotate('{:s}, {:5.1f} km'.format(k3, altx),
                                     [0, 1.09], xycoords='axes fraction',
                                     horizontalalignment='center',
                                     verticalalignment='center')
                        # colorbar
                        hcb = plt.colorbar(hc)#, ticks=np.arange(-50, 51, 10))
                        plt.tight_layout()
                        fig.savefig('/home/guod/WD4T/work3/'+
                                '{:s}_{:s}_{:d}_'.format(kk, k3, k2)+k0+
                                '_den_vortex.pdf')
                        plt.close(fig)
    return


def func6():
    # Find the minimum density and its location as a function of altitude (plot).
    # Density change as a function of neutral vorticity(scatter)
    # polar contours of density (fill) and vorticity (line)
    from spacepy.datamodel import dmarray
    import gitm
    import gitm_vorticity as gv
    from apexpy import Apex
    from matplotlib.ticker import AutoMinorLocator
    import gitm_3D_const_alt as g3ca
    timedate = pd.Timestamp('2010-03-23 06:00:00')
    if True:
        nlat, slat = 90, 60
    else:
        nlat, slat = -60, -90
    # Figure 1: minimum density versus altitude
    fig1 = plt.figure()
    ax1 = plt.subplot()
    # Figure 4: MLAT of density minimum
    fig4 = plt.figure()
    ax4 = plt.subplot()
    # Figure 5: MLT of density minimum
    fig5 = plt.figure()
    ax5 = plt.subplot()
    # Figure 2: density change versus vorticity
    salt = [200, 400]
    fig2, ax2 = plt.subplots(
            len(salt), 2, sharex=True, sharey='row',
            figsize=[len(salt)*3, len(salt)*3.2])
    # Figure 3: density and vorticity changes versus LT and lat
    fig3, ax3 = plt.subplots(
            len(salt), 2,  figsize=[len(salt)*3,
            len(salt)*3.2], subplot_kw={'polar':True})
    for k11, k1 in enumerate(['1', '2']):
        fn = ('/home/guod/WD4T/gitm/run_imfby/run' +
              k1+'b/data/3DALL_t'+timedate.strftime('%y%m%d_%H%M%S')+'.bin')
        g = gitm.GitmBin(fn, varlist=['Rho', 'V!Dn!N (north)',
                                      'V!Dn!N (east)', 'V!Dn!N (up)'])
        gv.calc_vorticity(g, 'neutral')
        g['vt'] = dmarray(
                np.sqrt(g['V!Dn!N (north)']**2+g['V!Dn!N (east)']**2),
                attrs={'units':'m/s', 'scale':'linear',
                       'name':'total wind velocity'})
        dlat, dlon, rho, vorticity, vt = (
                g[k] for k in ['dLat', 'dLon', 'Rho', 'nvorticity', 'vt'])
        # latitude and longitude limits
        ilat = (dlat[0, :, 0] >= slat) & (dlat[0, :, 0]<=nlat)
        ilon = (dlon[:, 0, 0] >= 0) & (dlon[:, 0, 0]<=360)
        # fetch dlat, dlon, rho, vorticity
        dlat, dlon, rho, vorticity, vt =(
                k[ilon][:, ilat] for k in (dlat, dlon, rho, vorticity, vt))
        rhomin = pd.DataFrame(columns=('alt', 'lat', 'lon', 'rhomin'))
        rhomax = pd.DataFrame(columns=('alt', 'lat', 'lon', 'rhomax'))
        for k0 in range(g.attrs['nAlt']):
            rhom = (np.mean(rho[..., k0]*np.cos(dlat[...,k0]/180*np.pi))
                    /np.mean(np.cos(dlat[..., k0]/180*np.pi)))
            rho[..., k0] = rho[..., k0]/rhom
            altt = g['Altitude'][0, 0, k0]
            # Minimum
            ind = np.argmin(rho[:, :, k0].reshape(-1))
            latt, lont, rhomint, vortt = (
                    k[..., k0].reshape(-1)[ind]
                    for k in [dlat, dlon, rho, vorticity])
            rhomint = pd.DataFrame(
                    np.array([altt, latt, lont, rhomint]).reshape(1, -1),
                    columns=('alt', 'lat', 'lon', 'rhomin'))
            rhomin = rhomin.append(rhomint)
            # Maximum
            ind = np.argmax(rho[:, :, k0].reshape(-1))
            latt, lont, rhomaxt, vortt = (
                    k[..., k0].reshape(-1)[ind]
                    for k in [dlat, dlon, rho, vorticity])
            rhomaxt = pd.DataFrame(
                    np.array([altt, latt, lont, rhomaxt]).reshape(1, -1),
                    columns=('alt', 'lat', 'lon', 'rhomax'))
            rhomax = rhomax.append(rhomaxt)
        apex=Apex(date=2010)
        mlat, mlon = apex.convert(
                rhomin.lat, rhomin.lon, source='geo', dest='qd',
                height=rhomin.alt/1000)
        mlatt, mlt = apex.convert(
                rhomin.lat, rhomin.lon, source='geo', dest='mlt',
                height=rhomin.alt/1000, datetime=timedate)
        # Figure 1, 4, 5
        ax1.plot(rhomin.rhomin, rhomin.alt/1000)
        ax4.plot(mlt, rhomin.alt/1000)
        ax5.plot(mlat, rhomin.alt/1000)
        ax1.set_xlim(0.5, 1)
        ax1.set_ylim(100, 700)
        ax1.grid()
        ax1.set_xlabel(r'$\rho_{min}/\rho_{mean}$', fontsize='x-large')
        ax1.set_ylabel('Altitude (km)')
        ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(10))
        ax1.legend(['Run1', 'Run2'])
        ax4.set_xlim(0, 12)
        ax4.set_xticks(np.arange(0, 13, 2))
        ax4.set_ylim(100, 700)
        ax4.grid()
        ax4.set_xlabel('MLT (hour)')
        ax4.set_ylabel('Altitude (km)')
        ax4.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax4.yaxis.set_minor_locator(AutoMinorLocator(10))
        ax4.legend(['Run1', 'Run2'])
        ax5.set_xlim(-90, 90)
        ax5.set_xticks(np.arange(-90, 91, 30))
        ax5.set_ylim(100, 700)
        ax5.grid()
        ax5.set_xlabel('MLAT (degree)')
        ax5.set_ylabel('Altitude (km)')
        ax5.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax5.yaxis.set_minor_locator(AutoMinorLocator(10))
        ax5.legend(['Run1', 'Run2'])

        # Figure 2
        dlat, dlon, rho, vorticity, vt = (
                g[k] for k in ['dLat', 'dLon', 'Rho', 'nvorticity', 'vt'])
        # latitude and longitude limits
        sslat = -85 if slat==-90 else slat
        nnlat = 85 if nlat==90 else nlat
        ilat = (dlat[0, :, 0] >= sslat) & (dlat[0, :, 0]<=nnlat)
        ilon = (dlon[:, 0, 0] >= 0) & (dlon[:, 0, 0]<=360)
        # fetch dlat, dlon, rho, vorticity
        dlat, dlon, rho, vorticity, vt =(
                k[ilon][:, ilat] for k in (dlat, dlon, rho, vorticity, vt))
        for k22, k2 in enumerate(salt):
            plt.sca(ax2[k22, k11])
            alt = g['Altitude'][0, 0, :]
            ialt = np.argmin(abs(alt-k2*1000)) # in GITM, unit of alt is m
            altx = alt[ialt]/1000
            rhot, vortt, vtt = (k[..., ialt] for k in (rho, vorticity, vt))
            rhot, vortt, vtt = (k.reshape(-1) for k in (rhot, vortt, vtt))
            ivtt = vtt>200
            rhot, vortt, vtt = (k[ivtt] for k in (rhot, vortt, vtt))
            plt.plot(vortt*1000, rhot,'k.')
        plt.xlim(-2, 2)
        [ax2[-1, k].set_xlabel('Vorticity (10$^{-3}$/s)') for k in range(2)]
        ax2[0, 0].set_title('Run1')
        ax2[0, 1].set_title('Run2')
        [ax2[k, 0].set_ylabel(r'$\rho$ (kg/m$^3$)') for k in range(len(salt))]
        plt.tight_layout()

        # Figure 3
        #rhorange=[[]]
        for k22, k2 in enumerate(salt):
            axt, hc = g3ca.contour_single(
                    ax3[k22, k11], 'Rho', 'polar', g, alt=k2, nlat=nlat,
                    slat=slat, nzlevels=20)
            #axt, hv =g3ca.vector_single(
            #        ax3[k22, k11], g, 'neu','polar', alt=k2, nlat=nlat,
            #        slat=slat)
            g3ca.contour_single(
                    ax3[k22, k11], 'nvorticity', 'polar', g,
                    alt=k2, nlat=nlat, slat=slat, zmin=-0.001, zmax=0.001,
                    nzlevels=11, fill=False,
                    colors='k', linewidths=1, alpha=0.5)
            #plt.colorbar(hc, ax=axt)
        ax3[0, 0].set_title('Run1', y=1.1)
        ax3[0, 1].set_title('Run2', y=1.1)
        #plt.tight_layout()


def func7():
    # Time constant of the density response to IMF By reversal

    NS = 'N'
    print('Results for '+NS+'H')
    fig, ax = plt.subplots(7, 1, sharex=True, figsize=(5.3, 8))
    # Time index only for imf by
    btime = pd.Timestamp('2010-03-20 00:00:00')
    etime = pd.Timestamp('2010-03-24 00:00:00')
    dt = pd.Timedelta('1min')
    time = pd.date_range(btime, etime, freq=dt)
    # Figure 1: IMF By
    plt.sca(ax[0])
    fn1 = '/home/guod/WD4T/gitm/run_imfby/run1c/imf1.dat'
    imf1 = pd.read_csv(
            fn1, delim_whitespace=True, comment='#',
            header=None,
            names=('year', 'month', 'day', 'hour', 'minute', 'second', 'ms',
                   'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
            usecols=['by'])
    imf1 = pd.DataFrame(np.array(imf1), index=time, columns=['By'])
    plt.plot(imf1.index, imf1.By)
    fn2 = '/home/guod/WD4T/gitm/run_imfby/run2c/imf2.dat'
    imf2 = pd.read_csv(
            fn2, delim_whitespace=True, comment='#',
            header=None,
            names=('year', 'month', 'day', 'hour', 'minute', 'second', 'ms',
                   'bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't'),
            usecols=['by'])
    imf2 = pd.DataFrame(np.array(imf2), index=time, columns=['By'])
    plt.plot(imf2.index, imf2.By)
    plt.ylabel(r'$B_y$')
    plt.ylim([-10, 10])

    # Find file names (100322_230000 -> 100323_060000)
    bdir = '/home/guod/big/raid4/guod/run_imfby/run1c/data/'
    tname = os.listdir(bdir)
    tname = [i for i in tname if '.bin' in os.path.splitext(i)]
    tname = np.sort(tname)
    inds, = np.where(tname=='3DALL_t100322_230000.bin')[0]
    #inds, = np.where(tname=='3DALL_t100322_000000.bin')[0]
    basename = tname[inds:]

    pdir = '/home/guod/big/raid4/guod/run_imfby/run2c/data/'
    tname = os.listdir(pdir)
    tname = [i for i in tname if '.bin' in os.path.splitext(i)]
    tname = np.sort(tname)
    inds, = np.where(tname=='3DALL_t100322_230000.bin')[0]
    pertname = tname[inds:]

    # Figure 2-5: ion, neutral wind vorticity and neutral density
    # Prepare data and save
    if False:
        print('Data preparation:')
        oo = pd.DataFrame()
        if True: # Run1
            bpdir = bdir
            bpname = basename
            savename = '/home/guod/data/tmp/w3_07_01.dat'
            print('Data preparation for Run 1')
        else: # Run2
            bpdir = pdir
            bpname = pertname
            savename = '/home/guod/data/tmp/w3_07_02.dat'
            print('Data preparation for Run 2')
        # limit latitudes
        nlat, slat = (90, 60) if NS=='N' else (-60, -90)
        print('Latitude range <{:d}, {:d}>'.format(slat, nlat))
        for k00, k0 in enumerate(bpname):
            # Read file
            gg = gitm.GitmBin(
                    bpdir+k0,
                    varlist=['Rho', 'V!Dn!N (north)', 'V!Dn!N (east)',
                             'V!Dn!N (up)', 'V!Di!N (north)', 'V!Di!N (east)',
                             'V!Di!N (up)', 'Temperature'])
            # Add vorticity
            gv.calc_vorticity(gg, neuion='neutral', name='nvorticity')
            gv.calc_vorticity(gg, neuion='ion', name='ivorticity')
            ilat = (gg['dLat'][0, :, 0]>slat) & (gg['dLat'][0, :, 0]<nlat)
            ilon = (gg['dLon'][:, 0, 0]>0) & (gg['dLon'][:, 0, 0]<360)
            for k1 in range(gg.attrs['nAlt']):
                alt, lat, lon, vorti, vortn, rho, Temp = (
                        gg[k][..., k1][:, ilat][ilon, ...].reshape(-1) \
                        for k in ['Altitude', 'dLat', 'dLon', 'ivorticity',
                                  'nvorticity', 'Rho', 'Temperature'])
                # find minimum neu density
                rhomini = np.argmin(rho)
                lat1, lon1 = lat[rhomini], lon[rhomini]
                #ilat1 = (lat < lat1+2.5) & (lat>lat1-2.5)
                ilat1 = (lat < lat1+0.01) & (lat>lat1-0.01)
                #ilon1 = (lon < lon1+5) & (lon>lon1-5)
                ilon1 = (lon < lon1+0.01) & (lon>lon1-0.01)
                rhomean = np.mean(rho*np.cos(lat/180*np.pi))\
                          / np.mean(np.cos(lat/180*np.pi))
                tempmean = np.mean(Temp*np.cos(lat/180*np.pi))\
                           / np.mean(np.cos(lat/180*np.pi))
                rhomin = np.mean(rho[(ilat1) & (ilon1)])
                rhot = rhomin/rhomean
                vortit = np.mean(vorti[(ilat1) & (ilon1)])
                vortnt = np.mean(vortn[(ilat1) & (ilon1)])
                oo = oo.append(
                        pd.DataFrame(
                            np.array([lat1, lon1, vortit, vortnt, rhot,
                                      rhomean, rhomin, tempmean]
                                    ).reshape(1, -1),
                            index=pd.MultiIndex.from_product(
                                [[gg['time']], [gg['Altitude'][0, 0, k1]]],
                                names=['time', 'alt']),
                            columns=['lat', 'lon', 'vorti', 'vortn',
                                     'rho', 'rhomean', 'rhomin', 'tempmean']))
        pd.to_pickle(oo, savename)
        print('End of data preparation.')
    fn = ('/home/guod/data/tmp/w3_07_01.dat',
          '/home/guod/data/tmp/w3_07_02.dat')
    alt = 400
    print('Altitude {:d} Km'.format(alt))
    for k in fn:
        from apexpy import Apex
        oo = pd.read_pickle(k)
        apex = Apex(date=2010)
        oo['mlat'], mlon = apex.convert(
                oo['lat'], oo['lon'], source='geo', dest='qd')
        mlat, oo['mlt'] = apex.convert(
                oo['lat'], oo['lon'], source='geo', dest='mlt',
                datetime=oo.index.get_level_values(0))
        altarray = oo.index.get_level_values(1).values
        ialt = np.argmin(np.abs(altarray-alt*1000))
        #oot = oo.loc[(slice(None), altarray[ialt]), :]
        oot = oo.xs(oo.index.get_level_values(1)[ialt], level=1)
        # ivorticity
        plt.sca(ax[1])
        vorti_rolling = oot['vorti'].rolling(
                window=1, min_periods=1, center=True).mean()
        plt.plot(oot.index, vorti_rolling*1e4)
        plt.ylabel('$V_i$')
        plt.ylim([-10, 20])
        # nvorticity
        plt.sca(ax[2])
        plt.plot(oot.index, oot.vortn*1e4)
        plt.ylabel('$V_n$')
        plt.ylim([0, 15])
        # minimum density
        plt.sca(ax[3])
        plt.plot(oot.index, oot.rhomin)
        plt.ylabel(r'$\rho_{min}$')
        # minimum density/mean density
        plt.sca(ax[4])
        plt.plot(oot.index, oot.rho)
        plt.ylabel(r'$\rho_{min}/\rho_{mean}$')
        plt.ylim([0.5, 1])
        # mean density
        plt.sca(ax[5])
        plt.plot(oot.index, oot.rhomean)
        plt.ylabel(r'$\rho_{mean}$')
        # mean Temperature
        plt.sca(ax[6])
        plt.plot(oot.index, oot.tempmean)
        plt.ylabel(r'$Temp_{mean}$')
        #    # mlat
        #    plt.sca(ax[6])
        #    plt.plot(oot.index, oot.mlat)
        #    plt.ylabel(r'MLAT')
        #    plt.ylim([60, 90])
        #    # mlat
        #    plt.sca(ax[7])
        #    plt.plot(oot.index, oot.mlt)
        #    plt.ylabel(r'MLT')
        #    plt.ylim([4, 10])
    # set axis
    plt.xlim('2010-03-22 23:00:00', '2010-03-23 06:00:00')
    plt.xticks(
            pd.date_range(
                '2010-03-22 23:00:00', '2010-03-23 06:00:00', freq='1H'),
            (pd.date_range(
                '2010-03-22 23:00:00', '2010-03-23 06:00:00', freq='1H')
             -pd.Timestamp('2010-03-23 00:00:00'))/pd.Timedelta('1H'))
    plt.xlabel(r'Hours from B$_y$ reversal')
    plt.tight_layout(h_pad=0.01)
    return


def func8():
    # Every 30 minutes, see the results.
    import gitm
    import gitm_3D_const_alt as g3ca
    from apexpy import Apex
    import gitm_divergence as gd
    stime = pd.Timestamp('2010-03-23 00:00:00')
    etime = pd.Timestamp('2010-03-23 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='30min')
    fn1 = ['/home/guod/big/raid4/guod/run_imfby/run1c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]
    fn2 = ['/home/guod/big/raid4/guod/run_imfby/run2c/data/3DALL_t'+
          k.strftime('%y%m%d_%H%M%S')+'.bin' for k in timeidx]
    #oo = pd.read_pickle('/home/guod/data/tmp/w3_05_02.dat')
    alt = 400
    path =  '/home/guod/Documents/work/fig/tmp/w03/'
    #altarray = oo.index.get_level_values(1).values
    #ialt = np.argmin(np.abs(altarray-alt*1000))
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)
    for k0 in range(len(fn1)):
        g1 = gitm.GitmBin(
                fn1[k0],
                varlist=['Rho', 'Temperature',
                         'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)'])
        g2 = gitm.GitmBin(
                fn2[k0],
                varlist=['Rho', 'Temperature',
                         'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)'])
        gd.calc_divergence(g2)
        # MLT of geomagnetic pole
        qlt = g2['time'].hour+g2['time'].minute/60+g2['time'].second/3600+qlon/15
        qtheta, qr = qlt/12*np.pi, 90-qlat
        # location of density minima
        #oot = oo.loc[(g2['time'], oo.index.get_level_values(1)[ialt]), :]
        #lat, lon = oot['lat'], oot['lon']
        #lt = g2['time'].hour+g2['time'].minute/60+g2['time'].second/3600+lon/15
        #theta, r = lt/12*np.pi, 90-lat

        # Density and wind run2
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'Rho', 'pol', g2, alt=alt, nlat=90, slat=60,
                zmin=4e-12, zmax=9e-12, nzlevels=30)
        ax, hv = g3ca.vector_single(
                ax, g2, 'neu', 'pol', alt=alt, nlat=90, slat=60,
                scale=1500, alpha=0.7)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_den_win_run2_'+g1['time'].strftime('%H%M')+'.pdf')

        # Density and wind run1
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'Rho', 'pol', g1, alt=alt, nlat=90, slat=60,
                zmin=4e-12, zmax=9e-12, nzlevels=30)
        ax, hv = g3ca.vector_single(
                ax, g1, 'neu', 'pol', alt=alt, nlat=90, slat=60,
                scale=1500, alpha=0.7)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_den_win_run1_'+g1['time'].strftime('%H%M')+'.pdf')

        # Density and wind difference
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_diff(
                ax, 'Rho', 'pol', g1, g2, alt=alt, nlat=90, slat=60,
                zmin=-30, zmax=30, nzlevels=30)
        ax, hv = g3ca.vector_diff(
                ax, g1, g2, 'neu', 'pol', alt=alt, nlat=90, slat=60,
                scale=1000, alpha=0.7)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_den_win_diff_'+g1['time'].strftime('%H%M')+'.pdf')

        # Temperature difference
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_diff(
                ax, 'Temperature', 'pol', g1, g2, alt=alt, nlat=90, slat=60,
                diff_type='absolute', nzlevels=30, zmin=-200, zmax=200)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_T_diff_'+g1['time'].strftime('%H%M')+'.pdf')

        # Temperature run1
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'Temperature', 'pol', g1, alt=alt, nlat=90, slat=60,
                nzlevels=20, zmin=1350, zmax=1650)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_T_run1_'+g1['time'].strftime('%H%M')+'.pdf')
        # Temperature run2
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'Temperature', 'pol', g2, alt=alt, nlat=90, slat=60,
                nzlevels=20, zmin=1350, zmax=1650)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_T_run2_'+g1['time'].strftime('%H%M')+'.pdf')
        # vertical wind difference
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_diff(
                ax, 'V!Dn!N (up)', 'pol', g1, g2, alt=alt, nlat=90, slat=60,
                diff_type='absolute', nzlevels=30, zmin=-30, zmax=30)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_vup_diff_'+g1['time'].strftime('%H%M')+'.pdf')
        # vertical wind run2
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'V!Dn!N (up)', 'pol', g2, alt=alt, nlat=90, slat=60,
                nzlevels=30, zmin=-50, zmax=50)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_vup_run2_'+g1['time'].strftime('%H%M')+'.pdf')
        # vertical wind run1
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'V!Dn!N (up)', 'pol', g1, alt=alt, nlat=90, slat=60,
                nzlevels=30, zmin=-50, zmax=50)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_vup_run1_'+g1['time'].strftime('%H%M')+'.pdf')
        # divergence run 2
        plt.figure()
        ax = plt.subplot(polar=True)
        ax, hc = g3ca.contour_single(
                ax, 'divergence', 'pol', g2, alt=alt, nlat=90, slat=60,
                nzlevels=20, zmin=-0.0005, zmax=0.0005)
        plt.colorbar(hc)
        #ax.scatter(theta, r, color='green')
        ax.scatter(qtheta, qr, color='k')
        plt.title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT')
        plt.savefig(path+'w03_func08_div_run2_'+g1['time'].strftime('%H%M')+'.pdf')
        plt.close('all')


def func9():
    import gitm
    # UT change in density
    dt = pd.date_range(start='2010-03-20 00:00:00', end='2010-03-23 06:00:00',
                       freq='1h')
    path = '/home/guod/big/raid4/guod/run_imfby/run1c/data/'
    fn = (path+'3DALL_t{:s}.bin'.format(k.strftime('%y%m%d_%H%M%S')) for k in dt)
    nlat, slat, alt = 90, 60, 400
    if False: # data preparation
        oo = pd.DataFrame(columns=('time', 'rhomean', 'tempmean'))
        for k00, k0 in enumerate(fn):
            if not os.path.isfile(k0):
                continue
            g = gitm.GitmBin(k0, varlist=['Rho', 'Temperature'])
            ilat = (g['dLat'][0, :, 0]>slat) & (g['dLat'][0, :, 0]<nlat)
            ilon = (g['dLon'][:, 0, 0]>=0) & (g['dLon'][:, 0, 0]<=360)
            ialt = np.argmin(np.abs(g['Altitude'][0, 0, :]-alt*1000))
            gt = g['time']
            rho, temp, lat = (g[k][ilon][:, ilat][:, :, ialt] for k in (
                              'Rho', 'Temperature', 'Latitude'))
            rhomean = np.mean(rho*np.cos(lat))/np.mean(np.cos(lat))
            tempmean = np.mean(temp*np.cos(lat))/np.mean(np.cos(lat))
            oo = oo.append(pd.DataFrame(
                [[gt, rhomean, tempmean]], columns=('time', 'rhomean', 'tempmean')))
        if path[-8] == '2':
            pd.to_pickle(oo, '/home/guod/data/tmp/w3_09_2.dat')
        if path[-8] == '1':
            pd.to_pickle(oo, '/home/guod/data/tmp/w3_09_1.dat')
    oo2 = pd.read_pickle('/home/guod/data/tmp/w3_09_2.dat')
    oo1 = pd.read_pickle('/home/guod/data/tmp/w3_09_1.dat')
    fig, ax = plt.subplots(2, 1, sharex=True)
    plt.sca(ax[0])
    plt.plot(oo1['time'], oo1['rhomean'], label='run1')
    plt.plot(oo2['time'], oo2['rhomean'], label='run2')
    plt.legend()
    plt.ylabel(r'$\rho$')
    plt.ylim(4.5*1e-12, 6*1e-12)
    plt.sca(ax[1])
    plt.plot(oo1['time'], oo1['tempmean'])
    plt.plot(oo2['time'], oo2['tempmean'])
    plt.ylabel(r'T')
    plt.xlim('2010-03-21 00:00:00', '2010-03-23 06:00:00')
    #plt.ylim(1300, 1500)
    xticks = pd.date_range('2010-03-21 00:00:00', '2010-03-23 06:00:00', freq='6H')
    plt.xticks(xticks, xticks.strftime('%H'))
    plt.xlabel('UT')
    return


# END
#------------------------------------------------------------
if __name__ == '__main__':
    plt.close('all')
    a = func7()
    plt.show()
