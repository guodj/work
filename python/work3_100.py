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
def polar_contourf_champ_grace(self, ax, whichcolumn='rho', ns='N', **kwargs):
    """ A contourf of multiple-day density versus time and latitude.

    Args:
        ax: axis handle
        whichcolumn: string, 'rho400', 'rho', 'rho410'.
        ns: northern or hemisphere
        **kwargs: for contourf
    Return:
        hc: handle of the contourf plot
    ----------------------------------------
    x axis: days from '2000-1-1'
    """
    from matplotlib.ticker import AutoMinorLocator
    if not self.empty:
        self['epochday'] = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
        btime = self['epochday'].min()
        etime = self['epochday'].max()

        self = self.add_updown()
        tmp = self[self.isup] if updown is 'up' else self[self.isdown]

        ut0 = np.arange(np.floor(btime), np.floor(etime)+1+0.1/24, 0.5/24)
        lat0 = np.arange(-90,91,3)
        ut, lat = np.meshgrid(ut0, lat0)
        rho = griddata((tmp['epochday'], tmp.lat),
                       tmp[whichcolumn], (ut, lat),
                       method='linear', rescale=True)
        for index, k in enumerate(ut0):
            fp = abs(tmp['epochday']-k)<0.5/24
            if not fp.any():
                rho[:,index]=np.nan

        hc = ax.contourf(ut, lat, rho, 10, **kwargs)

        ax.set_xlim(np.floor(btime),np.floor(etime)+1)
        ax.set_xticks(np.arange(np.floor(btime),np.floor(etime)+2))
        ax.set_xticklabels(pd.date_range(
                tmp.index[0],
                tmp.index[-1]+pd.Timedelta('1d')).
                strftime('%j'))
        ax.set_ylim(-90,90)
        ax.set_yticks(np.arange(-90,91,30))
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(3))
        ax.tick_params(which='both', width=1.2)
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)
        ax.set_title('LT: {:.1f}'.format(tmp['LT'].median()))
        ax.set_xlabel('Day of {:d}'
                      .format(tmp.index[0].year),fontsize=14)
        ax.set_ylabel('Latitude', fontsize=14)
        return hc#, rho
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
    hours = mdates.HourLocator(range(0,25,3))
    hoursfmt = mdates.DateFormatter('%H')
    fig1,ax1 = plt.subplots(3,1,sharex=True,figsize=(7,7)) # for IMF and AE
    plot_omni(ax1,'2009-11-1','2010-12-31',['Bym','Bzm','AE'],'5m')
    for k0 in range(2):
        plt.sca(ax1[k0])
        plt.ylim(-10,10)
        plt.yticks(np.arange(-10,11,5))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.grid(dashes=(4,1))
        plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
    plt.sca(ax1[2])
    plt.ylim(0,800)
    plt.yticks(np.arange(0,801,200))
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.grid(dashes=(4,1))
    ax1[-1].xaxis.set_major_locator(hours)
    ax1[-1].xaxis.set_major_formatter(hoursfmt)
    # For satellites
    if False:
        denchamp = get_champ_grace_data('2009-11-1','2010-12-31',satellite='champ')
        dengrace = get_champ_grace_data('2009-11-1','2010-12-31',satellite='grace')
        dengoce = get_goce_data('2009-11-1','2010-12-31')
        pd.to_pickle((denchamp, dengrace, dengoce),DATADIR+'tmp/w3_00.dat')
    denchamp, dengrace, dengoce = pd.read_pickle(DATADIR+'tmp/w3_00.dat')
    for date in pd.date_range('2009-11-1','2010-12-31'):
        ax1[-1].set_xlim(date,date+pd.Timedelta('2D'))
        ax1[-1].set_xlabel('Hours of date: '+
                           date.date().strftime('%Y-%m-%d')+
                           '/'+(date+pd.Timedelta('1D')).date().strftime('%d'))
        plt.tight_layout()
        # By default, figures won't change until end of the script
        # draw() forces a figure redraw
        draw()
        # for satellites
        denchamptmp = ChampDensity(denchamp[date:date+pd.Timedelta('2D')])
        dengracetmp = ChampDensity(dengrace[date:date+pd.Timedelta('2D')])
        dengocetmp = GoceData(dengoce[date:date+pd.Timedelta('2D')])
        fig2 = plt.figure(figsize=(7,6.5))
        ax2 = plt.subplot(polar=True)
        hc = denchamptmp.satellite_position_lt_lat(ns='N')
        hc.set_color('r') if not hc is None else None
        hc.set_label('CHAMP') if not hc is None else None
        hc = dengracetmp.satellite_position_lt_lat(ns='N')
        hc.set_color('b') if not hc is None else None
        hc.set_label('GRACE') if not hc is None else None
        hc = dengocetmp.satellite_position_lt_lat(ns='N')
        hc.set_color('k') if not hc is None else None
        hc.set_label('GOCE') if not hc is None else None
        #--------------------------------------------------------------------------------
        # Set polar(lat, LT) coordinates
        ax2.set_rmax(30)
        ax2.set_rgrids(np.arange(10,31,10),['$80^\circ$','$70^\circ$','$60^\circ$'],fontsize=14)
        ax2.set_theta_zero_location('S')
        ax2.set_thetagrids(np.arange(0,361,90),[0,6,12,18],fontsize=14,frac=1.05)
        ax2.set_title('Date: '+date.date().strftime('%Y-%m-%d')+'/'+
                      (date+pd.Timedelta('1D')).date().strftime('%d'))
        #--------------------------------------------------------------------------------
        hl = ax2.legend(loc=(0.85,0.8))
        plt.show()
        input()
        plt.close(fig2)
    return

def func2():
    #----------------------------------------
    # Check one of the interesting cases: 2010-5-29
    #----------------------------------------
    bdate = '2010-5-29'
    edate = '2010-5-29 21:0:0'
    mdate = '2010-5-29 12:0:0'
    dench = get_champ_grace_data(bdate,edate,satellite='champ')
    dench['marglat'] = mf.lat2arglat(dench.lat)
    dengr = get_champ_grace_data(bdate,edate,satellite='grace')
    dengr['marglat'] = mf.lat2arglat(dengr.lat)
    dengo = get_goce_data(bdate,edate)
    dengo['marglat'] = mf.lat2arglat(dengo.lat)
    den = (dengo,dench,dengr)
    # Calculate orbit number
    for k0 in range(3):
        diffarglat = np.insert(np.diff(den[k0].marglat), 0, 0)
        den[k0]['orbitn'] = 0
        den[k0].loc[diffarglat<0, 'orbitn']=1
        den[k0]['orbitn'] = den[k0]['orbitn'].cumsum()
    fig,ax = plt.subplots(3, 1, sharex=True, figsize=(5.4, 7))
    yl = ((-1.5e-11, 1.5e-11), (-1e-11, 1e-11), (-3e-13, 3e-13))
    boundinglat = 50
    for k0 in range(3):
        plt.sca(ax[k0])
        tmp = den[k0]
        mn = tmp.orbitn.max()
        for k1 in range(mn+1):
            tmp1 = tmp[tmp.orbitn == k1]
            tmp1 = tmp1[tmp1.lat>boundinglat] # Mid and high latitudes.
            # Below uses tmp1.index.min(), because thermosphere needs
            # time to respond to forcings.
            c = 'r' if tmp1.index.min()<pd.Timestamp(mdate) else 'b'
            plt.plot(tmp1.marglat, tmp1.rho-tmp1.rho.mean(),
                    linestyle='--', color=c, alpha=0.5)
        plt.ylim(yl[k0])
        plt.xlim(boundinglat,180-boundinglat)
    plt.show()
    return den

# END
if __name__ == '__main__':
    plt.close('all')
    a = func2()
