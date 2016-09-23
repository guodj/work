#--------------------------------------------------------------------------------
# For the 3rd project
#
# By Dongjie, USTC/UM, on Tue Sep 20 02:39:02 CST 2016
#--------------------------------------------------------------------------------

#Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from champ_grace import *
from goce import *
from omni import *

def func1():
    #----------------------------------------
    # Find interesting cases in 2010
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
    plot_omni(ax1,'2009-11-1','2010-12-31',['Bym','Bzm','AE'],'5minute')
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
        denchamp = get_champ_grace_data(pd.date_range('2009-11-1','2010-12-31'),satellite='champ')
        dengrace = get_champ_grace_data(pd.date_range('2009-11-1','2010-12-31'),satellite='grace')
        dengoce = get_goce_data(pd.date_range('2009-11-1','2010-12-31'))
        pd.to_pickle((denchamp, dengrace, dengoce),'/data/tmp/w3_00.dat')
    denchamp, dengrace, dengoce = pd.read_pickle('/data/tmp/w3_00.dat')
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
        dengocetmp = GoceDensity(dengoce[date:date+pd.Timedelta('2D')])
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

# END
if __name__ == '__main__':
    func1()
