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


def func4():
    import gitm
    import pandas as pd
    import gitm_3D_lat_lt as gll
    from apexpy import Apex
    fig, ax = plt.subplots(4, 2, figsize=[6, 13],
                           subplot_kw=dict(projection='polar'))
    champ = ChampDensity('2010-02-26', '2010-02-27')
    nlat, slat = -30, -90
    plot_type='polar'
    alt = 400 # km
    zkey = 'Rho'
    for k in range(3, 25, 3):
        plt.sca(ax[int((k/3-1)//2), int(k/3%2-1)])
        tt0 = pd.Timestamp('2010-02-26')+pd.Timedelta(k-3, 'h')
        tt1 = pd.Timestamp('2010-02-26')+pd.Timedelta(k, 'h')
        # add champ orbit
        champt = champ[tt0:tt1]
        champt.__class__ = ChampDensity
        champt.satellite_position_lt_lat(
                plt.gca(), nlat=nlat, slat=slat, plot_type=plot_type,
                zorder=100, color='gray', alpha=0.5)
        # add msis difference density
        ttstr0 = tt0.strftime('%y%m%d_%H%M%S')
        ttstr1 = tt1.strftime('%y%m%d_%H%M%S')
        a0 = gitm.GitmBin(
                '/home/guod/tmp/3DALL_msis/3DALL_t'+ttstr0+'.bin',
                varlist=[zkey])
        a1 = gitm.GitmBin(
                '/home/guod/tmp/3DALL_msis/3DALL_t'+ttstr1+'.bin',
                varlist=[zkey])
        axt, hcont = gll.gitm_3D_lat_lt_diff(
                plt.gca(), zkey, plot_type, a0, a1, alt=alt,
                diff_type='relative', title=False, figname=None, draw=True,
                nlat=nlat, slat=slat, dlat=30, dlt=6, lt00='S', zmax=15,
                zmin=-15, zcolor=None, data_type="contour")
        plt.title(tt1.strftime('%H')+'-'+tt0.strftime('%H'))
        # add magnetic pole
        csign =1 if nlat>0 else -1
        pp = csign*90
        mplat, mplon = Apex().convert(pp, 0, source='qd', dest='geo')
        mput = tt1.hour+tt1.minute/60+tt1.second/3600
        mplt = mput+mplon/15
        mpr, mptheta = 90-csign*mplat, mplt/12*np.pi
        plt.plot(mptheta, mpr, 'ko')
        # It is strange that the minimum radius is not 0
        axt.set_rmin(0)
    [ax[-1, k].set_xlabel('LT') for k in range(2)]
    plt.subplots_adjust(top=0.95, bottom=0.07)
    cax = plt.axes([0.3, 0.025, 0.4, 0.01])
    plt.colorbar(hcont, cax=cax, ticks=np.arange(-15, 16, 5),
                 orientation='horizontal')
    plt.show()

# END
#------------------------------------------------------------
if __name__ == '__main__':
    plt.close('all')
    a = func4()
    plt.show()
