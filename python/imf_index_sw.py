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
def get_imf(dates,dsource = 'Omin',resolution = '1hour'):
    """Obtain IMF data in selected dates

    Args:
        dates: pd.DatetimeIndex or array of strings.
        dsource: 'ACE' or 'Omin'
        resolution: '5minute' or '1hour'
    returns:
        pd.DataFrame, IMF data indexed with datetime
    """
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    if dsource == 'ACE':
        # ACE only has 1h resolution data
        fpath = '/data/ACE_IMF_1h/'
        fname = ['{:s}imf{:s}_csv.txt'.format(fpath,k) for k in years]
        bad_data = -999.9
        imf = [pd.read_csv(fn,parse_dates = [0],index_col = [0])
               for fn in fname if os.path.isfile(fn)]
    elif dsource == 'Omin':
        if resolution == '5minute':
            fpath = '/data/Omin_SW_IMF_5min/IMF_5min/csv/'
            fname = ['{:s}imf_5min_{:s}_csv.lst'.format(fpath,k) for k in years]
            bad_data = 9999.99
        elif resolution == '1hour':
            fpath = '/data/Omin_SW_IMF_1h/IMF/csv/'
            fname = ['{:s}imf{:s}_csv.txt'.format(fpath,k) for k in years]
            bad_data = 999.9
        imf = [pd.read_csv(fn,parse_dates = [0],index_col = [0])
               for fn in fname if os.path.isfile(fn)]
    if imf:
        imf = pd.concat(imf)
        fp = np.floor(imf.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
        imf = imf.loc[fp]
        imf = imf.replace(bad_data,np.nan)
        return imf
    else:
        return pd.DataFrame()


def get_index(dates):
    """Obtain solar or geomagnetic indics data in selected dates.
    Only 1-hour resolution data is included

    Args:
        dates: pd.DatetimeIndex or array of strings.
    returns:
        pd.DataFrame, data indexed with datetime
    """
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    fname = ['/data/Omin_Solar_Geo_index_1h/csv/'
             'index{:s}_csv.txt'.format(y) for y in years]
    index = [pd.read_csv(
            fn,
            parse_dates=[0],
            index_col=[0]) for fn in fname if os.path.isfile(fn)]
    if index:
        index = pd.concat(index)
        fp = np.floor(index.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
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
    fname = ['/data/Omin_SW_IMF_1h/SW_l1/csv/'
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


def plot_imf_index_sw(bdate, edate, variables):
    # Plot variables in panels
    dates = pd.date_range(bdate,edate)
    imf = get_imf(dates)
    index = get_index(dates)
    sw = get_sw(dates)
    comb = pd.concat([imf,index,sw],axis=1)

    nv = len(variables)
    fig, ax = plt.subplots(nv,1,sharex=True)
    for k00, k0 in enumerate(variables):
        plt.sca(ax[k00])
        plt.plot(comb.index, comb[k0])
        plt.title(k0)
    return fig, ax


#TEST
#--------------------------------------------------------------------------------
if __name__ == '__main__':
    from pylab import *
    import matplotlib.dates as mdates
    hours = mdates.HourLocator(range(0,25,3))
    hoursfmt = mdates.DateFormatter('%H')
    fig,ax = plot_imf_index_sw('2010-1-1','2010-12-31',['Bx','Bye'])
    for k0 in range(2):
        plt.sca(ax[k0])
        ax[k0].set_ylim(-10,10)
        plt.grid(axis='y',dashes=(4,1))
    ax[1].xaxis.set_major_locator(hours)
    ax[1].xaxis.set_major_formatter(hoursfmt)
    plt.show()
    for date in pd.date_range('2010-1-1','2010-12-31'):
        ax[1].set_xlim(date,date+pd.Timedelta('1D'))
        ax[1].set_xlabel(date.date())
        # By default, figures won't change until end of the script
        # draw() forces a figure redraw
        draw()
        input()
