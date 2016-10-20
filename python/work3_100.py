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
    #----------------------------------------
    # Find interesting cases in 2009 (11,12) and 2010
    # interesting case:
    #     1, IMF By from - to + or the opposite
    #     2, CHAMP, GRACE and GOCE orbits overlap
    #----------------------------------------
    # For IMF and AE
    from pylab import draw # Use draw()
    import matplotlib.dates as mdates
    from matplotlib.ticker import AutoMinorLocator
    if False:
        omni = get_omni('2009-11-1','2010-12-31', ['Bym', 'Bzm', 'V'], '5m')
        omni['EF'] = omni.V*(np.sqrt(omni.Bym*omni.Bym + omni.Bzm*omni.Bzm)-omni.Bzm)/2/1000
        pd.to_pickle(omni,DATADIR + 'tmp/w3_01.dat')
    omni = pd.read_pickle(DATADIR + 'tmp/w3_01.dat')
    hours = mdates.HourLocator(range(0,25,3))
    hoursfmt = mdates.DateFormatter('%H')
    fig1,ax1 = plt.subplots(3,1,sharex=True,figsize=(5,6.7)) # By, Bz, E
    # Case 1
    #    bdate = '2010-5-29'
    #    edate = '2010-5-30 00:00:00'
    # Case 2
    #    bdate = '2009-11-14 00:00:00'
    #    edate = '2009-11-15 00:00:00'
    # Case 3
    #    bdate = '2009-12-12 15:00:00'
    #    edate = '2009-12-13 06:00:00'
    # Case 4
    #    bdate = '2010-01-02 00:00:00'
    #    edate = '2010-01-03 00:00:00'
    # Case 5
    #    bdate = '2010-01-03 00:00:00'
    #    edate = '2010-01-04 00:00:00'
    # Case 6
    #    bdate = '2010-01-05 09:00:00'
    #    edate = '2010-01-06 03:00:00'
    # Case 7
    #    bdate = '2010-04-12 09:00:00'
    #    edate = '2010-04-13 00:00:00'
    # Case 8
    bdate = '2010-05-18 21:00:00'
    edate = '2010-05-20 00:00:00'
    pv = ['Bym', 'Bzm', 'EF']
    ylb = ['$B_y$ (nT)', '$B_z$ (nT)', 'Electric Field (mV/m)']
    yl = [(-15, 15), (-15, 15), (0, 10)]
    yt = [np.arange(-15, 16, 5), np.arange(-15, 16, 5), np.arange(0, 11, 2)]
    for k0 in range(3):
        plt.sca(ax1[k0])
        tmp = omni[pv[k0]][bdate:edate]
        plt.plot(tmp.index, tmp)
        plt.xlim(bdate, edate)
        plt.ylabel(ylb[k0])
        plt.ylim(yl[k0])
        plt.yticks(yt[k0])
        plt.grid(dashes=(4,1))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
    for k0 in range(2):
        plt.sca(ax1[k0])
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
    plt.sca(ax1[2])
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
    ax1[-1].xaxis.set_major_locator(hours)
    ax1[-1].xaxis.set_major_formatter(hoursfmt)
    plt.xlabel('Hours of {:s}'.format(pd.Timestamp(bdate).strftime('%Y-%m-%d')))
    plt.tight_layout()
    #    # For satellites
    #    if False:
    #        denchamp = get_champ_grace_data('2009-11-1','2010-12-31',satellite='champ')
    #        dengrace = get_champ_grace_data('2009-11-1','2010-12-31',satellite='grace')
    #        dengoce = get_goce_data('2009-11-1','2010-12-31')
    #        pd.to_pickle((denchamp, dengrace, dengoce),DATADIR+'tmp/w3_00.dat')
    #    denchamp, dengrace, dengoce = pd.read_pickle(DATADIR+'tmp/w3_00.dat')
    #    for date in pd.date_range('2009-11-1','2010-12-31'):
    #        ax1[-1].set_xlim(date,date+pd.Timedelta('2D'))
    #        ax1[-1].set_xlabel('Hours of date: '+
    #                           date.date().strftime('%Y-%m-%d')+
    #                           '/'+(date+pd.Timedelta('1D')).date().strftime('%d'))
    #        plt.tight_layout()
    #        # By default, figures won't change until end of the script
    #        # draw() forces a figure redraw
    #        draw()
    #        # for satellites
    #        denchamptmp = ChampDensity(denchamp[date:date+pd.Timedelta('2D')])
    #        dengracetmp = ChampDensity(dengrace[date:date+pd.Timedelta('2D')])
    #        dengocetmp = GoceData(dengoce[date:date+pd.Timedelta('2D')])
    #        fig2 = plt.figure(figsize=(7,6.5))
    #        ax2 = plt.subplot(polar=True)
    #        hc = denchamptmp.satellite_position_lt_lat(ns='N')
    #        hc.set_color('r') if not hc is None else None
    #        hc.set_label('CHAMP') if not hc is None else None
    #        hc = dengracetmp.satellite_position_lt_lat(ns='N')
    #        hc.set_color('b') if not hc is None else None
    #        hc.set_label('GRACE') if not hc is None else None
    #        hc = dengocetmp.satellite_position_lt_lat(ns='N')
    #        hc.set_color('k') if not hc is None else None
    #        hc.set_label('GOCE') if not hc is None else None
    #        #--------------------------------------------------------------------------------
    #        # Set polar(lat, LT) coordinates
    #        ax2.set_rmax(30)
    #        ax2.set_rgrids(np.arange(10,31,10),['$80^\circ$','$70^\circ$','$60^\circ$'],fontsize=14)
    #        ax2.set_theta_zero_location('S')
    #        ax2.set_thetagrids(np.arange(0,361,90),[0,6,12,18],fontsize=14,frac=1.05)
    #        ax2.set_title('Date: '+date.date().strftime('%Y-%m-%d')+'/'+
    #                      (date+pd.Timedelta('1D')).date().strftime('%d'))
    #        #--------------------------------------------------------------------------------
    #        hl = ax2.legend(loc=(0.85,0.8))
    #        plt.show()
    #        input()
    #        plt.close(fig2)
    #    return

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
    dench = get_champ_grace_density(bdate,edate,satellite='champ')
    dench['arglat'] = mf.lat2arglat(dench.lat)
    mdench = dench.groupby([np.floor(dench.arglat/3)*3, dench.index<mdate])['rho'].mean()
    dengr = get_champ_grace_density(bdate,edate,satellite='grace')
    dengr['arglat'] = mf.lat2arglat(dengr.lat)
    mdengr = dengr.groupby([np.floor(dengr.arglat/3)*3, dengr.index<mdate])['rho'].mean()
    dengo = get_goce_data(bdate,edate)
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
    plt.show()
    return mdench

# END
if __name__ == '__main__':
    plt.close('all')
    a = func1()
    plt.show()
