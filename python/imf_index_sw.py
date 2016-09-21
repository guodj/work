#--------------------------------------------------------------------------------
# Functions to handle IMF, index, solar wind parameters.
#
# By Dongjie, USTC, Sat Sep 17 09:17:16 CST 2016
#
# Contain:
#     get_imf: Get IMF data from Omin (5m or 1H) or ACE (1H),
#     get_index: Get Kp, ap, Dst, pc, AE, AU, AL, f107, R indices
#     get_sw: Get solar wind temperature, pressure, speed, proten
#             density, flow longitude, flow latitude
#     get_cirlist: Get CIR list created by McPherron
#     get_cirlist_noicme: Get CIR list with the file interfacenoicme
#     ......
#--------------------------------------------------------------------------------

# Global imports
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
def get_imf(dates,resolution = '1hour'):
    """Obtain Omin IMF data in selected dates

    Args:
        dates: pd.DatetimeIndex or array of strings.
        resolution: '5minute' or '1hour'
    returns:
        pd.DataFrame, IMF data indexed with datetime
    """
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    if resolution in ['5minute','5m']:
        fpath = '/data/Omin_5m/IMF_5min/csv/'
        fname = ['{:s}imf_5min_{:s}_csv.lst'.format(fpath,k) for k in years]
        bad_data = 9999.99
    elif resolution in ['1hour','1h']:
        fpath = '/data/Omin_1h/IMF/csv/'
        fname = ['{:s}imf{:s}_csv.txt'.format(fpath,k) for k in years]
        bad_data = 999.9
    imf = [pd.read_csv(fn,parse_dates = [0],index_col = [0])
           for fn in fname if os.path.isfile(fn)]
    if imf:
        imf = pd.concat(imf)
        # Two ways to convert date to number:
        # 1, to_julian_date; 2, (date1-date2)/'1D'
        fp = np.floor(imf.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
        imf = imf.loc[fp]
        imf = imf.replace(bad_data,np.nan)
        return imf
    else:
        return pd.DataFrame()


def get_index(dates,resolution = '1hour'):
    """Obtain solar or geomagnetic indics data in selected dates.

    Args:
        dates: pd.DatetimeIndex or array of strings.
    returns:
        pd.DataFrame, data indexed with datetime
    """
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    if resolution in ['1hour','1h','1H']:
        fname = ['/data/Omin_1h/index/csv/'
                 'index{:s}_csv.txt'.format(y) for y in years]
    elif resolution in ['5minute','5m']:
        # Only for AE and PC
        fname = ['/data/Omin_5m/IMF_AE_PC_01_12_5m/csv/'
                 '{:s}_csv.lst'.format(y) for y in years]
    index = [pd.read_csv(
            fn,
            parse_dates=[0],
            index_col=[0]) for fn in fname if os.path.isfile(fn)]
    if index:
        index = pd.concat(index)
        fp = np.floor(index.index.to_julian_date()+0.5).isin(
                dates.to_julian_date()+0.5)
        if resolution in ['1hour','1h','1H']:
            index = index.loc[fp]
            index.loc[index.Kp==99, 'Kp'] = np.nan
            index.loc[index.R==999, 'R'] = np.nan
            index.loc[index.Dst==99999, 'Dst'] = np.nan
            index.loc[index.ap==999, 'ap'] = np.nan
            index.loc[index.f107==999.9, 'f107'] = np.nan
            index.loc[index.AE==9999, 'AE'] = np.nan
            index.loc[index.AL==99999, 'AL'] = np.nan
            index.loc[index.AU==99999, 'AU'] = np.nan
            index.loc[index.pc==999.9, 'pc'] = np.nan
        elif resolution in ['5minute','5m']:
            index.rename(columns={'PC':'pc'},inplace=True)
            index.loc[index.AE==99999, 'AE'] = np.nan
            index.loc[index.pc==999.99, 'pc'] = np.nan
        return index
    else:
        return pd.DataFrame()


def get_sw(dates):
    """Obtain solar wind data in selected dates.
    Only 1-hour resolution data is included

    Args:
        dates: pd.DatetimeIndex or array of strings.
    returns:
        pd.DataFrame, data indexed with datetime
    """
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    fname = ['/data/Omin_1h/SW_l1/csv/'
             'plasma{:s}_csv.txt'.format(y) for y in years]
    plasma = [pd.read_csv(
            fn,
            parse_dates=[0],
            index_col=[0]) for fn in fname if os.path.isfile(fn)]
    if plasma:
        plasma = pd.concat(plasma)
        fp = np.floor(plasma.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
        plasma = plasma.loc[fp]
        plasma.loc[plasma.temperature==9999999., 'temperature'] = np.nan
        plasma.loc[plasma.proten_density==999.9, 'proten_density'] = np.nan
        plasma.loc[plasma.speed==9999., 'speed'] = np.nan
        plasma.loc[plasma.flow_lon==999.9, 'flow_lon'] = np.nan
        plasma.loc[plasma.flow_lat==999.9, 'flow_lat'] = np.nan
        plasma.loc[plasma.pressure==99.99, 'pressure'] = np.nan
        return plasma
    else:
        return pd.DataFrame()


def get_cirlist():
    """ Obtain the CIR list during 1995-2006 created by R. McPherron

    Return: pd.DatetimeIndex
    """
    fname = '/data/CIRlist/streaminterfacelist.txt'
    cirlist = pd.read_csv(
            fname, delim_whitespace=True, comment='%', header=None,
            usecols=[0,1], names=['date','time'],
            parse_dates={'datetime': [0,1]})
    return pd.DatetimeIndex(cirlist.datetime)


def get_cirlist_noicme():
    """ Obtain the CIR list during 1998-2009 with the file interfacenoicme

    Return: pd.DatetimeIndex
    """
    fname = '/data/CIRlist/interfacenoicme.txt'
    cirlist = pd.read_csv(
            fname, delim_whitespace=True, comment='%', header=None,
            usecols=[0,1], names=['date','time'],
            parse_dates={'datetime': [0,1]})
    return pd.DatetimeIndex(cirlist.datetime)


def plot_imf_index_sw(ax, bdate, edate, variables,res='1hour', **kwargs):
    # ax should be created with plt.subplots() outside
    # variables should be a list or tuple even if only one variable
    # kwargs are parameters for plot
    # Plot variables in panels
    # res can be '1h' or '5m', '5m' is only for IMF, AE, PC
    vl = {'Bx':'$B_x$ (nT)', 'Bye':'GSE $B_y$ (nT)', 'Bze':'GSE $B_z$ (nT)',
          'Bym':'GSM $B_y$ (nT)', 'Bzm':'GSM $B_z$ (nT)','AE':'AE','pc':'PC'}
    dates = pd.date_range(bdate,edate)
    imf = get_imf(dates,res)
    index = get_index(dates,res)
    sw = get_sw(dates)
    comb = pd.concat([imf,index,sw],axis=1)

    for k00, k0 in enumerate(variables):
        plt.sca(ax[k00]) if len(variables)>1 else plt.sca(ax)
        plt.plot(comb.index, comb[k0], **kwargs)
        plt.ylabel(vl[k0])
    return

#TEST
#--------------------------------------------------------------------------------
if __name__ == '__main__':
    from pylab import *
    import matplotlib.dates as mdates
    from matplotlib.ticker import AutoMinorLocator
    hours = mdates.HourLocator(range(0,25,3))
    hoursfmt = mdates.DateFormatter('%H')
    fig,ax = plt.subplots(1,1,sharex=True)
    plot_imf_index_sw(ax,'2010-1-1','2010-1-31',['Bx'],'5minute')
    #fig,ax = plt.subplots(3,1,sharex=True)
    #plot_imf_index_sw(ax,'2010-1-1','2010-12-31',['Bx','Bye','AE'],'5minute')
    #for k0 in range(2):
    #    plt.sca(ax[k0])
    #    plt.ylim(-10,10)
    #    plt.yticks(np.arange(-10,11,5))
    #    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
    #    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    #    plt.grid(dashes=(4,1))
    #    plt.axhline(0,color='r',linestyle='--',dashes=[4,1])
    #plt.sca(ax[2])
    #plt.ylim(0,800)
    #plt.yticks(np.arange(0,801,200))
    #plt.gca().xaxis.set_minor_locator(AutoMinorLocator(3))
    #plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
    #plt.grid(dashes=(4,1))
    #ax[-1].xaxis.set_major_locator(hours)
    #ax[-1].xaxis.set_major_formatter(hoursfmt)
    #plt.show()
    #for date in pd.date_range('2010-1-1','2010-12-31'):
    #    ax[-1].set_xlim(date,date+pd.Timedelta('1D'))
    #    ax[-1].set_xlabel('Hours of date: '+date.date().strftime('%Y-%m-%d'))
    #    # By default, figures won't change until end of the script
    #    # draw() forces a figure redraw
    #    draw()
    #    input()
