#-*- coding: utf-8 -*-

"""For the second work
"""

__author__ = 'Dongjie Guo'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import griddata
import pdb   # set breakpoint
import matplotlib.dates as mdates
from scipy.signal import argrelextrema
from apexpy import Apex
from mpl_toolkits.basemap import Basemap
import matplotlib
import time


def get_date_polarity():
    """ Get solar wind sector polarities and their dates.

    Args:
        no input
    Returns:
        date_polarity: pd.DataFrame with columns 'polarity', 'season'.
            The index is dates
    """

    sb_fname = '/data/SBlist/SBlist.txt'
    sblist = pd.read_csv(
        sb_fname,
        sep='\s+',
        parse_dates={'dates':['year', 'month', 'day']},
        index_col='dates')
    date_polarity = pd.DataFrame()
    alist = []
    for row in sblist.itertuples():
        # row = ('1926-01-29', '-,+', 21, 5)
        index = (row[0] + pd.TimedeltaIndex(range(-row[2], row[3]), 'd'))
        value = list(row[2]*row[1][0]+row[3]*row[1][2])
        df = pd.DataFrame(value, index=index, columns=['polarity'])
        alist.append(df)
    date_polarity = pd.concat(alist)
    date_polarity = date_polarity.groupby(level=0,sort=False).first()
    date_polarity.replace(['+','-'],['away','toward'],inplace=True)
    doy = date_polarity.index.dayofyear
    date_polarity.ix[(doy>35)  & (doy<125),'season'] = 'me'
    date_polarity.ix[(doy>221) & (doy<311),'season'] = 'se'
    date_polarity.ix[(doy>128) & (doy<218),'season'] = 'js'
    date_polarity.ix[(doy>311) | (doy<36),'season'] = 'ds'
    return date_polarity


def get_sblist():
    """ Get solar wind sector polarity reversing dates

    Args:
        No input
    Returns:
        A dataframe indexed by dates with columns sbtype, lday, rday and season.
    """
    sb_fname = '/data/SBlist/SBlist.txt'
    sblist = pd.read_csv(
        sb_fname,
        sep='\s+',
        parse_dates={'dates':['year','month','day']},
        index_col='dates')
    sblist.replace(['+,-','-,+'], ['away-toward','toward-away'], inplace=True)
    doy = sblist.index.dayofyear
    sblist.ix[(doy>=35)  & (doy<=125),'season'] = 'me'
    sblist.ix[(doy>=221) & (doy<=311),'season'] = 'se'
    sblist.ix[(doy>125) & (doy<221),'season'] = 'js'
    sblist.ix[(doy>311) | (doy<35),'season'] = 'ds'
    return sblist


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


def get_goce_data(dates):
    """ get goce data during specified dates.

    Args:
        dates: pd.DatetimeIndex or array of strings.
    Returns:
        dataframe of goce density indexed with datetime. the columns
        are: alt, long, lat, LT, arglat, rho, cr_wnd_e, cr_wnd_n, cr_wnd_u
    """
    dates = pd.to_datetime(dates)
    yearmonth = np.unique(dates.strftime('%Y-%m'))
    yearmonth =pd.to_datetime(yearmonth)
    fname = ['/data/GOCE/data/goce_denswind_v1_3_{:s}-{:s}.txt'.format(
        k.strftime('%Y'),
        k.strftime('%m')) for k in yearmonth]
    goce_data = [pd.read_csv(
        fn, delim_whitespace=True, header=None, comment='#',
        names=['date', 'time', 'Time_system',
               'alt', 'long', 'lat', 'LT', 'arglat', 'rho',
               'cr_wnd_e', 'cr_wnd_n', 'cr_wnd_u'],
        parse_dates={'datetime':[0,1]},index_col=0)
        for fn in fname if os.path.isfile(fn)]
    if goce_data:
        goce_data = pd.concat(goce_data)
        fp = np.floor(goce_data.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
        goce_data = goce_data.loc[fp]
        return goce_data
    else:
        return pd.DataFrame()

class ChampDensity(pd.DataFrame):
    """ Subclass of pd.DataFrame FOR champ or grace density data.
    """
    #----------------------------------------------------------------#
    # Another name for ChampDensity.lat
    def get_lat(self):
        return self.lat
    def set_lat(self, value):
        self.lat = value
    def del_lat(self):
        del self.lat
    latitude = property(get_lat, set_lat, del_lat,'Latitudes in the data')
    # Another name for ChampDensity.long
    def get_lon(self):
        return self.long
    def set_lon(self, value):
        self.long = value
    def del_lon(self):
        del self.long
    lon = property(get_lon, set_lon, del_lon,'longitude in the data')
    longitude = property(get_lon, set_lon, del_lon,'longitude in the data')
    # Another name for ChampDensity.height
    def get_height(self):
        return self.height
    def set_height(self, value):
        self.height = value
    def del_height(self):
        del self.height
    altitude = property(get_height, set_height, del_height,
                        'altitude in the data')
    alt = property(get_height, set_height, del_height,'altitude in the data')
    # Another name for ChampDensity.rho
    def get_rho(self):
        return self.rho
    def set_rho(self, value):
        self.rho = value
    def del_rho(self):
        del self.rho
    density = property(get_rho, set_rho, del_rho,'density in the data')
    # Another name for ChampDensity.rho400
    def get_rho400(self):
        return self.rho400
    def set_rho400(self, value):
        self.rho400 = value
    def del_rho400(self):
        del self.rho400
    density400 = property(get_rho400, set_rho400, del_rho400,
                          'density at 400km in the data')
    # Another name for ChampDensity.Mlong
    def get_Mlong(self):
        return self.Mlong
    def set_Mlong(self, value):
        self.Mlong = value
    def del_Mlong(self):
        del self.Mlong
    Mlon = property(get_Mlong, set_Mlong, del_Mlong,
                    'Magnetic longitude in the data')
    #-------------------------------------------------------------------------#
    def subset_season(self, season):
        doy = self.index.dayofyear
        if season in ['me', 'ME', 'Me']:
            self = self[(35<=doy) & (doy<=125)]
        if season in ['se', 'SE', 'Se']:
            self = self[(221<=doy) & (doy<=311)]
        if season in ['js', 'JS', 'Js']:
            self = self[(128<=doy) & (doy<=218)]
        if season in ['ds', 'DS', 'Ds']:
            self = self[(311<=doy) | (doy<=36)]
        return ChampDensity(self)

    def add_updown(self, whichlat='lat3'):
        """ add 'isup' and 'isdown' columns to self
        Note that the function is appropriate for continuous data

        Args:
            whichlat: data column name used to seperate up and down
                orbits
        Returns:
            self added with columns 'isup' and 'isdown'

        Note:
            The results may be wrong near poles.
            Some bugs may exist at the data gap.
        """
        if not self.empty:
            lat = self[whichlat]
            dlat = lat.diff()
            dlat.iloc[0] = dlat.iloc[1]
            if whichlat=='lat3':
                fp = (np.abs(dlat)!=3) & (np.abs(dlat)!=0)
                fpid = np.argwhere(fp).reshape(-1)
                dlat[fpid] = dlat[fpid+1]
            self['isup'] = (dlat > 0)
            self['isdown'] = (dlat < 0)
            if whichlat=='lat3':
                fp1 = self.lat3>=87
                fp2 = dlat==0
                fp = fp1 & fp2
                self.loc[fp,'isdown'] = True
                fp1 = self.lat3<=-87
                fp2 = dlat==0
                fp = fp1 & fp2
                self.loc[fp,'isup'] = True
            self = self[(self.isup) | (self.isdown)]
            return ChampDensity(self)


    def get_lt_lat_density(self, whichdensity='rho400', latwindow=10,
                           ltwindow=8, interpolate=False, latn=180, ltn=48):
        """ Get lt, lat, data for polar contourf. This function is for a large
        data set that covers all the lt and lat ranges.

        Args:
            whichdensity: rho, rho400, rho410...
            latwindow: lat window, harmonic of 90: 3,10,30....
            ltwindow: local time window, harmonic of 24: 1,2,3,4,6,8...
            interpolate: bool, whether to interpolate data to higher resolution
            latn, ltn: for interpolation, bin numbers of lat and lt, integral
                multiple of 180 and 24 respectively
        Returns:
            lt, lat, density. All are 2D array

        Note: the sample number in each lat*lt grid is not considered
              maybe the interpolation is unnecessary
        """
        df = self[['lat', 'LT', whichdensity]]
        df['lat'] = pd.cut(df['lat'], bins=np.arange(-90,90+1e-10,latwindow),
                     labels=np.arange(-90,90+1e-10,latwindow)[:-1]+latwindow/2,
                     include_lowest=True)
        df['LT'] = pd.cut(df['LT'], bins=np.arange(0,24+1e-10,ltwindow),
                     labels=np.arange(0,24+1e-10,ltwindow)[:-1]+ltwindow/2,
                     include_lowest=True)
        df = df.groupby(['lat', 'LT']).mean().reset_index().pivot(
            'lat', 'LT', whichdensity)
        df[df.columns[0]+24] = df[df.columns[0]] # data is periodic in lt
        if interpolate:
            lt = np.linspace(df.columns.min(), df.columns.max(), ltn)
            lat = np.linspace(-90, 90, latn)
            lt, lat = np.meshgrid(lt, lat)
            df = df.unstack()
            df.name = whichdensity
            df = df.reset_index()
            data = griddata((df['lat'], df['LT']/12*180),
                            df[ whichdensity], (lat, lt/12*180), method='cubic')
            return lt, lat, data
        else:
            lt, lat = np.meshgrid(df.columns, df.index)
            return lt, lat, df.values


    def get_lon_lat_density(self, whichdensity='rho400', latwindow=10,
                            lonwindow=45, interpolate=False,
                            latn=180, lonn=360):
        """ return lon, lat, density for polar contourf. This function is for a
        large data set that covers all the lon and lat ranges.

        Args:
            whichdensity: rho, rho400, rho410...
            latwindow: lat window, harmonic of 90: 3,10...
            lonwindow: lon window, harmonic of 360: 10,30...
            interpolate: bool, whether to interpolate data to higher resolution
            latn, lonn: for interpolation, bin numbers of  lat and lon,
                intergral multiple of 180 and 360, respectively.
        Returns:
            lon, lat, density

        Note: the sample number in each lat*lon grid is not considered
              maybe the interpolation is unnecessary
        """
        df = self[['lat', 'lon', whichdensity]]
        df['lat'] = pd.cut(
                df['lat'], bins=np.arange(-90,90+1e-10,latwindow),
                labels=np.arange(-90,90+1e-10,latwindow)[:-1]+latwindow/2,
                include_lowest=True)
        df['lon'] = pd.cut(
                df['lon'], bins=np.arange(-180,180+1e-10,lonwindow),
                labels=np.arange(-180,180+1e-10,lonwindow)[:-1]+lonwindow/2,
                include_lowest=True)
        df = df.groupby(['lat', 'lon']).mean().reset_index().pivot(
            'lat', 'lon', whichdensity)
        df[df.columns[0]+360] = df[df.columns[0]] # data is periodic in lon
        if interpolate:
            lon = np.linspace(df.columns.min(), df.columns.max(), lonn)
            lat = np.linspace(-90, 90, latn)
            lon, lat = np.meshgrid(lon, lat)
            df = df.unstack()
            df.name = whichdensity
            df = df.reset_index()
            data = griddata((df['lat'], df['lon']),
                            df[ whichdensity], (lat, lon), method='cubic')
            return lon, lat, data
        else:
            lon, lat = np.meshgrid(df.columns, df.index)
            return lon, lat, df.values


    def LT_median(self):
        """ Get the local time of the ascending and descending satellite orbit.
        Only for a short-period continuous period.

        Returns:
            upLT, downLT: Local times for ascending and descending orbits.
        Note:
            if local times are [0,0,24,24], the result is 12, this is wrong.
        """
        if not self.empty:
            self = self.add_updown()
            grouped = self.groupby(self.isup)['LT']
            lt = grouped.median()
            return lt[True], lt[False]
        return np.nan, np.nan


    def orbit_mean(self,lats=(-85,85),updown='up'):
        """ Get the orbit mean density during specified latitudes and
        specified ascending or descending orbit

        Input:
            lats: two elements tuple, selected latitudes: lats[0]<=lat3<=lats[1]
            updown: ascending or descending orbit

        Output:
            result: DataFrame, columns: longitude, height, LT, rho, rho400

        Note:
            longitudeand LT may be wrong if lats are high latitudes
        """
        self = self.add_updown()
        tmp = self[self.isup] if updown =='up' else self[~self.isup]  # ascending or descending orbit?
        tmp = tmp[(tmp.lat3>=lats[0]) & (tmp.lat3<=lats[1])]  #  which latitudes?
        tmp['float_time'] = (
                tmp.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
        # Below is a good method to calculate orbit number
        tmp['difft'] = np.insert(
                np.diff(tmp.index)/pd.Timedelta('1h'), 0, np.nan)
        tmp['orbitn'] = 0
        tmp.loc[tmp.difft>0.5,'orbitn']=1
        tmp['orbitn'] = tmp['orbitn'].cumsum()

        # There are at least half of the expected points, only for CHAMP and GRACE data
        grouped = tmp.groupby('orbitn').filter(
                lambda x: len(x.index)>(lats[1]-lats[0])/3/2).groupby('orbitn')
        result = grouped.agg(
                {'float_time': np.nanmean,
                 'long':np.nanmedian,
                 'height':np.nanmean,
                 'LT':np.nanmedian,
                 'rho':np.nanmean,
                 'rho400':np.nanmean})
        result['float_time'] = (result['float_time']*24*3600).round()
        result['datetime'] = (pd.Timestamp('2000-1-1') +
                              pd.TimedeltaIndex(result.float_time,'S'))
        result = result.set_index('datetime')
        result = result.drop('float_time',axis=1)
        return result


    def scatter_lt_lat(self, upmarker='o', downmarker='o'):
        """ Show the lt and lat positions of the satellite in a polar
        coordinate.

        Input:
            upmarker: the same as 'marker' parameter in plt.scatter
            downmarker: the same as 'marker' parameter in plt.scatter

        Output:
            hcup, hcdown: scatter handles for the up and down orbits,
                respectively.
        """
        if not self.empty:
            self = self.add_updown()
            thetaup = self.loc[self.isup, 'LT']/12*np.pi
            rup = 90 - self.loc[self.isup, 'lat']
            hcup = plt.scatter(thetaup, rup, marker=upmarker, linewidths=0)

            thetadown = self.loc[self.isdown, 'LT']/12*np.pi
            rdown = 90 - self.loc[self.isdown, 'lat']
            hcdown = plt.scatter(thetadown, rdown,
                                 marker=downmarker, linewidths=0)
            return hcup, hcdown


    def contourf_date_lat(self, ax, whichcolumn='rho400',
                          updown='up', **kwargs):
        """ A contourf of multiple-day density versus date and latitude.

        Args:
            ax: axis handle
            whichcolumn: string, 'rho400', 'rho', 'rho410'.
            updown: string, 'up' or 'down'
            **kwargs: for contourf
        Return:
            hc: handle of the contourf plot
        """
        if not self.empty:
            self['epochday'] = ((self.index-self.index.min())
                    .total_seconds()/(24*3600))
            btime = self['epochday'].min()
            etime = self['epochday'].max()

            self = self.add_updown()
            if updown == 'up':
                tmp = self[self.isup]
            elif updown == 'down':
                tmp = self[self.isdown]

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
                    strftime('%m-%d'),rotation=45)
            ax.set_ylim(-90,90)
            ax.set_yticks(np.arange(-90,91,30))
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_minor_locator(AutoMinorLocator(3))
            ax.tick_params(which='both', width=1.2)
            ax.tick_params(which='major', length=7)
            ax.tick_params(which='minor', length=4)
            ax.set_title('LT: {:.1f}'.format(tmp['LT'].median()))
            ax.set_xlabel('Date of Year: {:d}'
                          .format(tmp.index[0].year),fontsize=14)
            ax.set_ylabel('Latitude', fontsize=14)
            return hc#, rho


    def relative_density_to_previous_orbit_lt(self, whichcolumn='rho400'):
        """Get the time interval (dhour_po) between two continuous orbits and
        density relative to in the previous orbit (rrho_po).
        (rho-rho_po)/(rho_po*dhour), unit: /h

        Input:
            whichcolumn: 'rho', 'rho400', 'rho410' ...

        Output:
            self with columns 'dhour_po', 'drho_po',
            'rrho_po' and 'distance_lt' in LT-lat coordinates added

        Note:
            The pole values may be unreliable
        """
        self = self.add_updown()
        groupdensity = self.groupby(['isup','lat3'])
        def rtl(x):
            x['dhour_po'] = np.insert(np.diff(x.index),0,'nat')/np.timedelta64(1,'h')
            x['rrho_po'] = x[whichcolumn].diff()/np.roll(x[whichcolumn],1)/x['dhour_po']
            x['drho_po'] = x[whichcolumn].diff()/x['dhour_po']

            lat = x['lat']/180*np.pi
            latpo = np.roll(lat, 1)
            deltaLT = x['LT'].diff()/12*np.pi
            x['distance_lt'] = 6378*np.arccos(
                    np.sin(lat)*np.sin(latpo) +
                    np.cos(lat)*np.cos(latpo)*np.cos(deltaLT))
            return x
        self = groupdensity.apply(rtl)
        fp1 = (self.dhour_po>1.6)
        fp2 = (self.dhour_po<1.45)
        fp3 = (self.distance_lt>600)
        fp = fp1 | fp2 | fp3
        self.loc[fp,'dhour_po'] = np.nan
        self.loc[fp,'rrho_po'] = np.nan
        self.loc[fp,'drho_po'] = np.nan
        self.loc[fp,'distance_lt'] = np.nan
        return ChampDensity(self)


    def relative_density_to_previous_orbit_mlt(self, whichcolumn='rho400'):
        """Relative density between 2 adjecent orbits.

        Input:
            whichcolumn: 'rho', 'rho400', 'rho410' ...

        Output:
            self with columns 'dhour_po', 'distance_mlt', 'rrho_po' and 'dpoints' added

        Note:
            The pole values may be unreliable
        """
        Mlat = np.array(self.Mlat)/180*np.pi
        Mlt = np.array(self.MLT)/12*np.pi
        time = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1h')
        rho = np.array(self[whichcolumn])

        dib, dl = 110, 16
        dd = np.ones([len(Mlat),dl])*np.nan
        dh = np.ones([len(Mlat),dl])*np.nan
        dr = np.ones([len(Mlat),dl])*np.nan
        for di in np.arange(dib,dib+dl):
            Mlatpo, Mltpo = np.roll(Mlat,di), np.roll(Mlt,di)
            timepo, rhopo = np.roll(time,di), np.roll(rho,di)
            dd[di:,di-dib] = 6378*np.arccos(
                    np.sin(Mlat[di:])*np.sin(Mlatpo[di:]) +
                    np.cos(Mlat[di:])*np.cos(Mlatpo[di:])*np.cos(Mlt[di:]-Mltpo[di:]))
            dh[di:,di-dib] = time[di:]-timepo[di:]
            dr[di:,di-dib] = (rho[di:]-rhopo[di:])/rhopo[di:]/dh[di:,di-dib]
        ddmin = np.argmin(dd,axis=1)

        ddout = dd[np.arange(len(Mlat)),ddmin]
        dhout = dh[np.arange(len(Mlat)),ddmin]
        drout = dr[np.arange(len(Mlat)),ddmin]
        dpout = ddmin + dib
        self['dhour_po'], self['distance_mlt'], self['rrho_po'], self['dpoints'] = (
                dhout, ddout, drout, dpout)
        fp1 = self.distance_mlt>800
        fp2 = (self.dhour_po>2) | (self.dhour_po<1)
        fp = fp1 | fp2
        self.loc[fp,['dhour_po', 'distance_mlt', 'rrho_po', 'dpoints']] = np.nan
        return ChampDensity(self)


    def add_index(self):
        """Add solar or geomagnetic parameters to ChampDensity.

        Output:
            self with columns 'ap', 'f107'... added
        """
        dates = pd.date_range(self.index.min().date(),
                              self.index.max().date(),
                              freq='1d')
        index = get_index(dates)
        index = index.reindex(self.index, method='ffill')
        self = pd.concat([self, index], axis=1)
        return ChampDensity(self)


    def add_imf(self):
        """Add IMF components to ChampDensity.

        Output:
            self with columns 'Bx', 'Bye'... added
        """
        dates = pd.date_range(self.index.min().date(),
                              self.index.max().date(),
                              freq='1d')
        imf = get_imf(dates)
        imf = imf.reindex(self.index, method='ffill')
        self = pd.concat([self, imf], axis=1)
        return ChampDensity(self)


    def add_epoch_sb(self, exclude_cir=True):
        """Add the epoch days relative to the solar wind sector polarity
        reversing dates.

        Output:
            self with columns 'epochat', 'epochta' added
        """
        sblist = get_sblist()
        sblist = sblist[self.index.min().date():self.index.max().date()]
        if exclude_cir:
            cirlist = get_cirlist()
            sblist = sblist[cirlist.min().date():
                            cirlist.max().date()]
            #cirlist = (
            #        cirlist.
            #        append([cirlist+pd.Timedelta('1D'),
            #                cirlist+pd.Timedelta('-1D')]).unique())
            sblist = sblist[~sblist.index.isin(pd.to_datetime(cirlist.date))]
        groupsb = sblist.groupby('sbtype')
        for sbat in groupsb.get_group('away-toward').itertuples():
            ldate = sbat[0]-pd.Timedelta(sbat[2],'d')
            rdate = sbat[0]+pd.Timedelta(sbat[3],'d')
            self.loc[ldate:rdate,'epochat'] = (
                    self[ldate:rdate].index - sbat[0])/np.timedelta64(1,'D')
        for sbta in groupsb.get_group('toward-away').itertuples():
            ldate = sbta[0]-pd.Timedelta(sbta[2],'d')
            rdate = sbta[0]+pd.Timedelta(sbta[3],'d')
            self.loc[ldate:rdate,'epochta'] = (
                    self[ldate:rdate].index - sbta[0])/np.timedelta64(1,'D')
        return ChampDensity(self)


    def add_epochtime(self, epoch0, lday=5, rday=5):
        """ add epoch times to self
        *** be careful of the overlap ****

        Input:
            epoch0: pd.DatatimeIndex, zero epoch of special cases
            lday & rday: the maximum days before and after 0 epoch time
        Output:
            self with column epochtime (unit: day)  added
        """
        for k in epoch0:
            maxlepoch = k - pd.Timedelta(lday,'d')
            maxrepoch = k + pd.Timedelta(rday,'d')
            self.loc[maxlepoch:maxrepoch, 'epochtime'] = (
                    self[maxlepoch:maxrepoch].index - k)/np.timedelta64(1,'D')
        return ChampDensity(self)

#------------------------end class-------------------------------------#


def set_lt_lat_polar(ax, pole='N', boundinglat=60):
    """Change some default sets of a LT-latitude polar coordinates

    Input
    ax: axis handle
    pole: 'N' or 'S', Northern or Southern Hemisphere
    """
    ax.set_theta_zero_location('S')
    ax.set_thetagrids(angles=[0,90,180,270], labels=[0,6,12,18], fontsize=14)
    if pole =='N':
        ax.set_rgrids(radii=[30,60,90],
                     labels=['$60^\circ$', '$30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    if pole =='S':
        ax.set_rgrids(radii=[30,60,90],
                     labels=['$-60^\circ$', '$-30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    ax.set_rmax(90-np.abs(boundinglat))


def set_lon_lat_polar(ax, pole='N',boundinglat=60):
    """change some default sets of a longitude-latitude polar coordinates

    Input
    ax: axis handle
    pole: 'N' or 'S', Northern or Southern Hemisphere
    """
    ax.set_theta_zero_location('S')
    ax.set_thetagrids(
            (-90,0,90,180),
            labels=('$-90^\circ$','$0^\circ$', '$90^\circ$','$\pm180^\circ$'),
            fontsize=14)
    if pole =='N':
        ax.set_rgrids(radii=[30,60,90],
                     labels=['$60^\circ$', '$30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    if pole =='S':
        ax.set_rgrids(radii=[30,60,90],
                     labels=['$-60^\circ$', '$-30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    ax.set_rmax(90-np.abs(boundinglat))


relative = lambda x: 100*(x-x.mean())/x.mean()


def get_density_dates(dates,satellite='champ'):
    """ get champ or grace data during specified dates.

    Args:
        dates: specified dates which can be args of pd.to_datetime(). and if a
            single date is given , it shoule be in format such as ['2005-1-1'],
            so pd.to_datetime() can convert it to DatetimeIndex class
        satellite: 'champ' or 'grace'
    Returns:
        dataframe of champ or grace density indexed with datetime. the columns
        are:lat3, lat, long, height, LT, Mlat, Mlong, MLT, rho, rho400,rho410,
        msis_rho, ...
    """
    dates = pd.to_datetime(dates)
    dates = dates[(dates>'2001-1-1') & (dates<'2011-1-1')]
    # champ and grace data are in this date range
    if satellite == 'champ':
        fname = ['/data/CHAMP23/csv/{:s}/ascii/'
                 'Density_3deg_{:s}_{:s}.ascii'.format(
                     k.strftime('%Y'),
                     k.strftime('%y'),
                     k.strftime('%j')) for k in dates]
    elif satellite == 'grace':
        fname = ['/data/Grace23/csv/{:s}/ascii/'
                 'Density_graceA_3deg_{:s}_{:s}.ascii'.format(
                     k.strftime('%Y'),
                     k.strftime('%y'),
                     k.strftime('%j')) for k in dates]
    rho = [pd.read_csv(
        fn,
        parse_dates=[0],
        index_col=[0]) for fn in fname if os.path.isfile(fn)]
    if rho:
        rho = pd.concat(rho)
        # Exclude duplicate points
        rho = rho.groupby(rho.index).first()  # pd.DataFrame.drop_duplicates() has something wrong
        return ChampDensity(rho)
    else:
        return ChampDensity()


def lt_lat_contourf(ax, lt, lat, data, whichhemisphere='N', **kwargs):
    """ Draw a polar contourf with lt as theta, lat as r

    Args:
        ax: axis handle
        lt, lat, data: local time, latitude and data, all are 2D arrays
        whichhemisphere: 'N' or 'S', Northern or Southern Hemisphere
        **kwargs: contourf args
    Return:
        hc: handle of the contourf plot
    """
    theta = lt/12*np.pi
    r = 90-lat if whichhemisphere=='N' else 90+lat
    hc = plt.contourf(theta, r, data, **kwargs)
    return hc


def lon_lat_contourf(ax, lon, lat, data, whichhemisphere='N', **kwargs):
    """ Draw a polar contourf with longitude as theta, latitude as r

    Args:
        ax: axis handle
        lon, lat, data: longitude, latitude and data, all are 2D arrays
        whichhemisphere: 'N' or 'S', Northern or Southern Hemisphere
        **kwargs: contourf args
    Return:
        hc: handle of the contourf plot
    """
    theta = lon/180*np.pi
    r = 90-lat if whichhemisphere=='N' else 90+lat
    hc = ax.contourf(theta, r, data, **kwargs)
    return hc


def get_satellite_lt():
    """ Get champ and grace LTs in every day during 2001-2010
    """
    import gc
    dates = pd.date_range('2001-1-1','2010-12-31')
    datelt = pd.DataFrame(columns=[ 'LTup', 'LTdown'])
    datelt = datelt.reindex(dates)
    for kdate in datelt.index:
        density = get_density_dates([kdate])
        datelt.loc[kdate,['LTup', 'LTdown']] = density.LT_median()
    datelt.index.name = 'date'
    datelt.to_csv('/data/CHAMP23/LT.dat',na_rep=np.nan)

    datelt1 = pd.DataFrame(columns=[ 'LTup', 'LTdown'])
    datelt1 = datelt1.reindex(dates)
    for kdate in datelt1.index:
        density = get_density_dates([kdate],satellite='grace')
        datelt1.loc[kdate,['LTup', 'LTdown']] = density.LT_median()
    datelt1.index.name = 'date'
    datelt1.to_csv('/data/Grace23/LT.dat',na_rep=np.nan)
    return


def great_circle_distance(loc1,loc2,r=6378):
    """calculate the great circle distance of two points on a sphere

    Input:
        loc1: [lat, lon] or (lat, lon), in degree
        loc2: [lat, lon] or (lat, lon)
        r: radius of the sphere

    Output:
        d: distance
    """
    lat1 = loc1[0]/180*np.pi
    lon1 = loc1[1]/180*np.pi
    lat2 = loc2[0]/180*np.pi
    lon2 = loc2[1]/180*np.pi
    d = r*np.arccos(np.sin(lat1)*np.sin(lat2) +
                    np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    return d


if __name__=='__main__':
    def f3():
        """ Variations of the CHAMP and GRACE daily average densities during 2001-2010
        return: plot
        x: date(2001-2010)
        y: F107, CHAMP and GRACE data
        """
        fpath = ('/data/tmp/t1.dat')
        if False:
            density1 = get_density_dates(pd.date_range('2001-1-1','2010-12-31'),
                                        satellite='champ')
            density2 = get_density_dates(pd.date_range('2001-1-1','2010-12-31'),
                                        satellite='grace')
            density1 = density1.resample('1D').mean()
            density1 = density1[['rho','rho400','rho410']]
            density2 = density2.resample('1D').mean()
            density2 = density2[['rho','rho400','rho410']]
            pd.to_pickle([density1,density2], fpath)
            import gc
            gc.collect()
        density1, density2 = pd.read_pickle(fpath)
        fig, ax = plt.subplots()
        plt.plot(density1['rho400']/1e-11,'ro',markeredgewidth=0, label='CHAMP',
                 zorder=1)
        plt.plot(density2['rho400']/1e-11,'bo',markeredgewidth=0, label='GRACE',
                 zorder =2)
        plt.xlim([pd.to_datetime('2001-1-1'),pd.to_datetime('2011-1-1')])
        ax.xaxis.set_major_locator(mdates.YearLocator(1,1,1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        plt.legend(frameon=False,loc=(0.7,0.8))
        plt.xlabel('Year',fontsize=14)
        plt.ylabel('Density at 400 km($10^{-11} kg/m^3$)',fontsize=14)
        plt.grid()

        index = get_index(pd.date_range('2001-1-1','2010-12-31'))
        index = index.resample('1D').mean()
        import gc
        gc.collect()
        ax2 = ax.twinx()
        ax2.plot(index['f107'],'o',color='grey',markersize=3,linewidth=1,
                 markeredgewidth=0,zorder=3,alpha=0.8,label='F10.7')
        ax2.set_ylim([-100,300])
        ax2.set_yticklabels([])
        ax2.legend(loc=(0.7,0.74),frameon=False)
        plt.show()
        return

    def f5():
        """  density response to cir, test relative_density_to_previous_orbit_lt
        """
        satellite = 'champ'
        lepoch = 2
        repoch = 7
        deltat = 0.5/24   #  epoch time bin window

        fpath = ('/data/tmp/t2.dat')
        if False:  # data preparation
            cirlist = get_cirlist_noicme()
            diffcirlist = np.insert(
                    np.diff(cirlist),0,'NaT')/np.timedelta64(1,'D')
            fp = np.roll(diffcirlist,-1)>=repoch
            cirlist = cirlist[fp]
            cirlist = cirlist[(cirlist>='2008-1-1') & (cirlist<='2008-12-31')]
            density = []
            for k in cirlist:
                rho = get_density_dates(k+pd.TimedeltaIndex(np.arange(-lepoch, repoch),'D'))
                rho = rho.relative_density_to_previous_orbit_lt('rho400')
                rho['epochtime'] = (rho.index - k)/pd.Timedelta('1D')
                density.append(rho)
            density = pd.concat(density)
            density[['lat3','epochtime','rrho_po', 'drho_po']].to_pickle(fpath)
            import gc
            gc.collect()
        density = pd.read_pickle(fpath)
        density['epochtime'] = pd.cut(
                density['epochtime'],
                bins=np.arange(-lepoch,repoch+1e-10,deltat),
                include_lowest=True,
                labels=np.arange(-lepoch,repoch,deltat)+deltat/2).astype('float')
        grouped = density.groupby(['lat3','epochtime'])

        fig, ax = plt.subplots(2,1,sharex=True)
        plt.sca(ax[0])
        tmp = grouped['rrho_po'].median().reset_index().pivot(
                index='lat3',columns='epochtime',values='rrho_po')
        tmp = (tmp*deltat*24+1)

        tmp1 = tmp.cumprod(axis='columns',skipna=False)
        tmp2 = (tmp.loc[:,tmp.columns<=-1]).prod(axis='columns')
        tmp_cumprod = (tmp1.divide(tmp2,axis='index')-1)*100

        lat3 = tmp_cumprod.index
        epochtime = tmp_cumprod.columns
        epochtime, lat3 = np.meshgrid(epochtime, lat3)
        plt.contourf(epochtime, lat3, tmp_cumprod.values,
                     levels=np.linspace(0,75,11))
        plt.xlim(-lepoch-1,repoch+1)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,30))
        plt.ylabel('Lat (deg)',fontsize=14)
        plt.colorbar(label=r'$\Delta\rho\ (\%)$')

        plt.sca(ax[1])
        tmp = grouped['drho_po'].median().reset_index().pivot(
                index='lat3',columns='epochtime',values='drho_po')
        tmp = (tmp*deltat*24)

        tmp1 = tmp.cumsum(axis='columns',skipna=False)
        tmp2 = (tmp.loc[:,tmp.columns<=-1]).sum(axis='columns')
        tmp_cumsum = tmp1.sub(tmp2,axis='index')

        lat3 = tmp_cumsum.index
        epochtime = tmp_cumsum.columns
        epochtime, lat3 = np.meshgrid(epochtime, lat3)
        plt.contourf(epochtime, lat3, tmp_cumsum.values/1e-13,
                     levels=np.linspace(0,3,11))
        plt.xlim(-lepoch-1,repoch+1)
        plt.xlabel('Epoch days', fontsize=14)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,30))
        plt.ylabel('Lat (deg)',fontsize=14)
        plt.colorbar(label=r'$\Delta\rho\ (10^{-13} kg/m^{-3})$')

        plt.tight_layout()
        plt.show()


    def f6():
        """ another method for calculating density response to CIR
            the result shows that this method may be better
        """
        ldays = 3
        rdays = 7
        fpath = '/data/tmp/t3.dat'
        if False:  # data preparation
            cirlist = get_cirlist_noicme()
            diffcirlist = np.insert(
                    np.diff(cirlist),0,'NaT')/np.timedelta64(1,'D')
            fp = np.roll(diffcirlist,-1)>=rdays
            cirlist = cirlist[fp]
            cirlist = cirlist[(cirlist>'2008-1-1') & (cirlist<'2009-1-1')]
            density = []
            def f(x):
                xb = x.loc[x.epochcir<=-1].ix[-1,'rho400']
                x['rrho'] = 100*(x['rho400']-xb)/xb
                return x
            for datetime in cirlist:
                date =datetime.date()
                dates = pd.Timestamp(date) + pd.TimedeltaIndex(
                        np.arange(-ldays,rdays+2),'D')
                tmp = get_density_dates(dates)
                tmp['epochcir'] = (tmp.index - datetime)/pd.Timedelta('1D')
                tmp = ChampDensity(tmp[(tmp.epochcir>-ldays) &
                                       (tmp.epochcir<rdays)])

                tmp = tmp.add_updown()
                grouped = tmp.groupby(['lat3','isup'])
                tmp = grouped.apply(f)

                density.append(tmp)
            density = pd.concat(density)
            density.to_pickle(fpath)
        density = pd.read_pickle(fpath)
        density = density[['lat3','epochcir','rrho']]
        density['epochcir'] = np.floor(density['epochcir']*16)/16
        density_m = (density.groupby(['lat3','epochcir']).
                     median().
                     reset_index().
                     pivot(index='lat3',columns='epochcir',values='rrho'))
        lat3 = density_m.index
        epochcir = density_m.columns
        epochcir,lat3 = np.meshgrid(epochcir,lat3)

        fig =plt.figure()
        plt.contourf(epochcir,lat3,density_m.values,levels=np.linspace(0,75,11))
        plt.xlim(-ldays,rdays)
        plt.xlabel('Epoch days')
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,30))
        plt.ylabel('Lat (deg)',fontsize=14)
        plt.colorbar(label=r'$\Delta\rho\ (\%)$')
        plt.show()


    def f7():
        """ CIR distribution with SB epoch days
        """
        sblist = get_sblist()
        sblist = sblist['2001-5-15':'2010-7-1'] # CHAMP date range
        cirlist = get_cirlist_noicme()
        cirlist = pd.DatetimeIndex(cirlist.date)

        cirdis = []
        sbepoch = np.arange(-30,31)
        for k in sbepoch:
            cirtmp = cirlist + pd.Timedelta(k,'D')
            cirdis.append(len(set(cirtmp) & set(sblist.index)))
        plt.bar(-sbepoch-0.4,cirdis,linewidth=0,color='k')  # why 0.4? refer to plt.bar?
        plt.xlim(sbepoch.min()-0.5,sbepoch.max()+0.5)
        plt.xticks(np.arange(-27,28,9),fontsize=14)
        plt.yticks(range(0,81,20),fontsize=14)
        plt.grid()
        plt.xlabel('Epoch Day',fontsize=14)
        plt.ylabel('No. of CIRs',fontsize=14)
        #plt.plot(sbepoch,cirdis)
        plt.show()


    def f8():
        """ find SBs without CIR
        A SB is selected when no CIR occurs during SB epoch days(-2~3).
        """
        sblist = get_sblist()
        sblist = sblist['2001-5-15':'2010-7-1'] # CHAMP date range
        cirlist = get_cirlist()
        cirlist = pd.DatetimeIndex(cirlist.date)
        sblist = sblist[(sblist.index>=cirlist.min()) &
                        (sblist.index<=cirlist.max())]
        cirlist_ext = cirlist
        for k in np.arange(-2, 4):
            cirlist_ext = cirlist_ext.append(cirlist + pd.Timedelta(k,'D'))
        cirlist_ext = cirlist_ext.unique()
        sblist = sblist[~sblist.index.isin(cirlist_ext)]
        if False:  # Show the result.
            a = sblist[sblist.sbtype=='away-toward'].groupby('season').apply(len)
            plt.bar(np.arange(4)+0.1, a,linewidth=0,color='k')
            plt.xlim(0,4)
            plt.xticks(np.arange(0.5,4.5),a.index)
        return sblist


    def f9():
        """ Superposed epoch analysis
        density variations during 10 sector polarity reversal centered days
        return: contourf
        x: epoch day
        y: latitude
        c: relative density variation
        """
        fpath = '/data/tmp/t4.dat'
        if False:  # data preparation
            satellite = 'champ'
            rho = []
            for sbtype in ['away-toward', 'toward-away']:
                sblist = get_sblist()
                sblist = f8()  #  SBs that are not accompanied with CIRs
                sblist = sblist[(sblist.sbtype==sbtype)]
                lst =[]
                for date, season in zip(sblist.index, sblist.season):
                    dates = date+pd.TimedeltaIndex(range(-2,3),'d')
                    rho_10days = get_density_dates(dates,satellite)
                    if not rho_10days.empty:
                        rho_10days['season'] = season
                        rho_10days = rho_10days.add_updown()
                        rho_10days = rho_10days[rho_10days.lat3.abs() <
                                rho_10days.lat3.max()]
                        grouped = rho_10days.groupby(['lat3','isup'])
                        def tmpf(x):
                            x['rrho400'] = 100*((x['rho400']-x['rho400'].mean())/
                                    x['rho400'].mean())
                            return x
                        rho_10days = grouped.apply(tmpf)
                        rho_10days.index = rho_10days.index - date
                        rho_10days = rho_10days[['lat3','rrho400','season']]
                    lst.append(rho_10days)
                rho.append(pd.concat(lst))
            pd.to_pickle(rho,fpath)

        density = pd.read_pickle(fpath)
        fig,ax=plt.subplots(4,2,sharex=True,sharey=True,figsize=[11,7])
        for k1,sbtype in enumerate(['away-toward','toward-away']):
            rho = density[k1]
            #rho = rho[(rho.long>-90) & (rho.long<0)]
            for k2, season in enumerate(['me','se','js','ds']):
                rhotmp = rho[rho.season == season]
                str_season=['Mar.','Sep.','Jun.','Dec.']
                plt.sca(ax[k2,k1])
                rhotmp['days'] = rhotmp.index.total_seconds()/(3600*24)
                rhotmp['days'] = pd.cut(
                        rhotmp['days'],
                        np.arange(-5,6,3/24),
                        labels=np.arange(-5,6,3/24)[:-1]+1.5/24,
                        include_lowest = True).astype('float')
                """longitudinal distribution in a grid
                rhotmp = rhotmp[(rhotmp.lat>80) & (rhotmp.days == 0.75/24)]
                plt.hist(rhotmp.long, 12,  facecolor='k', alpha=0.75)
                plt.xlim(-180,180)
                plt.xticks(np.arange(-180,181,60))
                ax[k2,0].set_ylabel('data points')
                ax[-1,k1].set_xlabel('longitude')
                """
                """line plot for high latitude
                rhotmp = rhotmp[rhotmp.lat < -60]
                xy = rhotmp.groupby('days')['rrho400'].median()
                plt.plot(xy.index, xy.values)
                """
                xyz = rhotmp.groupby(['lat3','days'])['rrho400'].median()
                xyz = xyz.reset_index()
                xyz = xyz.pivot('lat3','days','rrho400')
                x0 = xyz.columns.values
                y0 = xyz.index.values
                z = xyz.values
                x,y = np.meshgrid(x0,y0)

                plt.contourf(x,y,z,levels=np.arange(-20,21,1), extend='both')
                plt.xlim(-2,3)
                plt.xticks(np.arange(-2,4))
                plt.ylim(-90,90)
                plt.yticks(np.arange(-90,91,30))
                plt.title(str_season[k2]+', '+sbtype)
                ax[k2,0].set_ylabel('Lat',fontsize=14)
                ax[-1,k1].set_xlabel('Epoch days',fontsize=14)
        fig.subplots_adjust(left=0.1,right=0.89,wspace=0.1)
        cax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
        cbar=plt.colorbar(cax=cax,ticks=np.arange(-20,21,4))
        cbar.set_label(r'$\rho$ (%)',fontsize=14)
        plt.show()
        return


    def f10():
        """Case study of the density response to SBs
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        nn = 0
        sdates = ['2002-9-15','2002-10-14','2003-9-30','2003-10-22',
                  '2006-10-7','2006-11-4','2007-2-22','2007-6-29',
                  '2008-10-18','2008-10-30','2009-9-13']
        sdates = pd.DatetimeIndex(sdates)
        for k in sdates:
            adates = k+pd.TimedeltaIndex(np.arange(-5,5),'D')
            sw = get_sw(adates)
            imf = get_imf(adates)
            sgindex = get_index(adates)
            rhoom = [1,1,1,1,1,1,1,1]
            ltom = [2,2,2,2,2,2,2,2]
            for k11,k1 in enumerate(['champ','grace']):
                rho = get_density_dates(adates,satellite=k1)
                if rho.empty:
                    break
                for k22,k2 in enumerate(['north','south']):
                    if k2 is 'north':
                        lats = (60,90)
                    if k2 is 'south':
                        lats = (-90,-60)
                    for k33,k3 in enumerate(['up', 'down']):
                        if k3 is 'up':
                            updown = 'up'
                        if k3 is 'down':
                            updown = 'down'
                        rhorho = rho.orbit_mean(lats=lats,updown=updown)
                        rhoom[k11*4+k22*2+k33] = rhorho
                        ltom[k11*4+k22*2+k33] = rhorho.LT.median()
            if rho.empty:
                continue
            else:
                nn = nn + 1
                print(nn)
            fig,ax = plt.subplots(7,1,sharex=True,figsize=(7,10))
            ax[0].plot(sw.index, sw.speed,'darkgreen',linewidth=2)
            ax[0].set_ylabel(r'$\mathsf{V_{SW}}$',fontsize=14)

            ax[1].plot(imf.index, imf.Bx,'darkblue',linewidth=2)
            ax[1].plot(imf.index, imf.Bye,'darkred',linewidth=2)
            ax[1].axhline(0,color='gray',linestyle='--',zorder=0)
            ax[1].set_ylabel(r'$B_x, B_y$',fontsize=14)

            ax[2].plot(imf.index, imf.Bzm,'darkcyan',linewidth=2)
            ax[2].axhline(0,color='gray',linestyle='--',zorder=0)
            ax[2].set_ylabel(r'$B_z$',fontsize=14)

            ax[3].plot(sgindex.index, sgindex.f107,'darkorange',linewidth=2)
            ax[3].set_ylabel(r'$F10.7$',fontsize=14)

            ax[4].plot(sgindex.index, sgindex.ap,'purple',linewidth=2)
            ax[4].set_ylabel(r'$ap$',fontsize=14)

            for k1,k2,k3 in zip([0,1,4,5],
                                ['b-','g-','r-','k-'],
                                ['CHAMP','CHAMP','GRACE','GRACE']):
                ax[5].plot(rhoom[k1].index, rhoom[k1].rho400/1e-12,k2,
                           linewidth=2,label=k3+'LT: {:.1f}'.format(ltom[k1]),zorder=1)
            plt.sca(ax[5])
            plt.legend(loc=0,frameon=False,fontsize=10,ncol=2)
            for k1,k2,k3 in zip([2,3,6,7],
                                ['b-','g-','r-','k-'],
                                ['CHAMP','CHAMP','GRACE','GRACE']):
                ax[6].plot(rhoom[k1].index, rhoom[k1].rho400/1e-12,k2,
                           linewidth=2,label=k3+'LT: {:.1f}'.format(ltom[k1]),zorder=1)
            plt.sca(ax[6])
            plt.legend(loc=0,frameon=False,fontsize=10,ncol=2)
            for k1 in range(7):
                ax[k1].set_xlim(adates.min(),adates.max()+pd.Timedelta('1D'))
                ax[k1].xaxis.set_ticks(adates)
                ax[k1].grid(axis='x')
                ax[k1].yaxis.set_major_locator(MaxNLocator(nbins=4,prune='upper'))
            ax[-1].xaxis.set_ticklabels(adates.strftime('%y/%m/%d'),rotation=45)
            plt.subplots_adjust(hspace=0.03,top=0.95,right=0.95)
            plt.show()
            input()
            plt.close()
        return


    def f11():
        """test griddata and contourf
        """
        density = get_density_dates(pd.date_range('2005-1-1','2005-1-10'))
        density.add_updown()
        density = density[density.isup]
        density = density[abs(density.lat3)<87]
        x = (density.index-pd.Timestamp('2004-12-31'))/pd.Timedelta('1D')
        y = density.lat
        z = density.rho400

        fig, ax = plt.subplots(3,2,sharex=False, sharey=True)
        plt.sca(ax[0,0])
        plt.tricontourf(x,y,z,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.title('plt.tricontourf')
        plt.ylabel('Lat')

        xg, yg = np.meshgrid(np.arange(0,10,0.021),np.arange(-90,90,0.1))
        plt.sca(ax[0,1])
        zg = griddata((x*10,y),z,(xg*10,yg),method='linear',rescale=False)
        plt.contourf(xg*100,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(100,1000)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.title('contourf x: x*10')
        # dense grids
        xg, yg = np.meshgrid(np.arange(0,10,0.021),np.arange(-90,90,0.1))
        plt.sca(ax[1,0])
        zg = griddata((x*10,y),z,(xg*10,yg),method='linear',rescale=False)
        plt.contourf(xg,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.ylabel('$\Delta$t=0.5H, $\Delta$lat=0.1')
        plt.title('griddata x: x*10')

        plt.sca(ax[1,1])
        zg = griddata((x/10,y),z,(xg/10,yg),method='linear',rescale=False)
        plt.contourf(xg,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.title('griddata x: x/10')

        xg, yg = np.meshgrid(np.arange(0,10),np.arange(-90,90,10))
        plt.sca(ax[2,0])
        zg = griddata((x*10,y),z,(xg*10,yg),method='linear',rescale=False)
        plt.contourf(xg,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.ylabel('$\Delta$t=1D, $\Delta$lat=10')
        plt.title('griddata x: x*10')
        plt.xlabel('Day of 2005')

        plt.sca(ax[2,1])
        zg = griddata((x/10,y),z,(xg/10,yg),method='linear',rescale=False)
        plt.contourf(xg,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.title('griddata x: x/10')
        plt.xlabel('Day of 2005')

        plt.tight_layout()
        plt.show()
        return

    def f12():
        """x axis: day of year; yaxis: epoch day
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']

        champlt = pd.read_csv('/data/CHAMP23/LT.dat',
                              index_col=[0],parse_dates=[0])
        gracelt = pd.read_csv('/data/Grace23/LT.dat',
                              index_col=[0],parse_dates=[0])
        if True: # test champlt and gracelt
            fig, ax = plt.subplots(2,1,sharex=True)
            plt.hold(True)
            ax[0].scatter(champlt.index, champlt.LTup, linewidth=0,color='b', label='up')
            ax[0].scatter(champlt.index, champlt.LTdown, linewidth=0,color='r', label='down')
            ax[0].set_title('CHAMP')
            ax[0].legend()
            ax[1].scatter(gracelt.index, gracelt.LTup, linewidth=0,color='b', label='up')
            ax[1].scatter(gracelt.index, gracelt.LTdown, linewidth=0,color='r', label='down')
            ax[1].set_title('GRACE')
            ax[1].set_xlabel('Year')
            for k in np.arange(2):
                ax[k].set_xlim(['2001-1-1','2011-1-1'])
                ax[k].set_ylim([0,24])
                ax[k].set_yticks(range(0,25,4))
            plt.show()
        # select specified sblist and LT
        groupsblist = sblist.groupby('sbtype')
        atlist = groupsblist.get_group('away-toward')

        fig = plt.figure()
        lrday = 5
        if False: #  data preparation
            a = []
            for k in ['champ','grace']:
                if k is 'champ':
                    tmp = champlt[champlt.index.isin(atlist.index)]
                else:
                    tmp = gracelt[gracelt.index.isin(atlist.index)]
                tmp1 = tmp[(tmp.LTup>12) | (tmp.LTup<12)]
                tmp2 = tmp[(tmp.LTdown>12) | (tmp.LTdown<12)]
                for k1 in ['up','down']:
                    tmp3 = tmp1 if k1 is 'up' else tmp2
                    for k2 in tmp3.index:
                        rho = get_density_dates(
                                k2+pd.TimedeltaIndex(np.arange(-lrday,lrday),'D'),
                                satellite=k)
                        rhot = rho.orbit_mean(lats=(-30,30),updown=k1)
                        #  exclude points with Kp>40
                        #smindex = get_index(k2+pd.TimedeltaIndex(np.arange(-lrday,lrday),'D'))
                        #smindex = smindex.reindex(rhot.index,method='ffill')
                        #rhot = rhot[smindex.Kp<=40]
                        #  select longitudes near the south pole
                        #fp, = argrelextrema(np.array(abs(rhot['long']+55)),np.less,order=5)
                        #rhot = rhot.iloc[fp]
                        print(rhot.shape[0])
                        if rhot.shape[0]>=128:
                            rhot['rrho'] = (100*(rhot['rho400']-rhot['rho400'].mean())/
                                    rhot['rho400'].mean())
                            rhot['doy'] = rhot.index.dayofyear
                            rhot['epochday'] = (rhot.index-k2)/pd.Timedelta('1D')
                            a.append(rhot)
            a = pd.concat(a)
            a.to_pickle('/data/tmp/t5.dat')
        a = pd.read_pickle('/data/tmp/t5.dat')
        x = a.doy
        y = a.epochday
        z = a.rrho
        x0, y0 = np.meshgrid(np.arange(1,366,10), np.arange(-5,5,1.5/24))
        z0 = griddata((x/5000,y),z,(x0/5000,y0))
        plt.contourf(x0,y0,z0,levels=np.linspace(-40,60,20))
        plt.colorbar()
        plt.show()
        return a


    def f13():
        """ champ and grace densities during epoch days of sblist
        output:
        data[dataframe,dataframe] in which data[0] for away-toward, data[1] for toward-away
        data1[dataframe,dataframe] for high-lats in the Northern hemisphere
        data2[dataframe,dataframe] for middle-lats in the Northern hemisphere
        data3[dataframe,dataframe] for low-lats
        data4[dataframe,dataframe] for middle-lats in the Southern hemisphere
        data5[dataframe,dataframe] for high-lats in the Southern hemisphere
        """
        import gc
        lrday = 5
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2011-1-1']
        data = [pd.DataFrame(),pd.DataFrame()]
        data1 = [pd.DataFrame(),pd.DataFrame()]
        data2 = [pd.DataFrame(),pd.DataFrame()]
        data3 = [pd.DataFrame(),pd.DataFrame()]
        data4 = [pd.DataFrame(),pd.DataFrame()]
        data5 = [pd.DataFrame(),pd.DataFrame()]
        nn = [0,0]
        for k1, k2 in enumerate(['away-toward', 'toward-away']):
            sbtmp = sblist[sblist.sbtype==k2]
            for idate in sbtmp.index:
                datelist = idate + pd.TimedeltaIndex(np.arange(2*lrday)-lrday, 'd')
                for satellite in ['champ','grace']:
                    density = get_density_dates(datelist,satellite=satellite)
                    if density.empty:
                        print('no data around', idate)
                        continue
                    if len(np.unique(density.index.dayofyear))!=2*lrday:
                        print('there is data gap around', idate)
                        continue
                    density = density.add_updown()
                    for updown in ['up','down']:
                        nn[k1] = nn[k1]+1
                        print(k2,': ', nn)
                        d = density.orbit_mean(updown=updown)
                        d['rrho400'] = 100*(d['rho400'] - d['rho400'].mean())/d['rho400'].mean()
                        d['epochday'] = (d.index - idate)/pd.Timedelta('1D')

                        d1 = density.orbit_mean(lats=(59,85),updown=updown)
                        d1['rrho400'] = 100*(d1['rho400'] - d1['rho400'].mean())/d1['rho400'].mean()
                        d1['epochday'] = (d1.index - idate)/pd.Timedelta('1D')

                        d2 = density.orbit_mean(lats=(29,61),updown=updown)
                        d2['rrho400'] = 100*(d2['rho400'] - d2['rho400'].mean())/d2['rho400'].mean()
                        d2['epochday'] = (d2.index - idate)/pd.Timedelta('1D')

                        d3 = density.orbit_mean(lats=(-31,31),updown=updown)
                        d3['rrho400'] = 100*(d3['rho400'] - d3['rho400'].mean())/d3['rho400'].mean()
                        d3['epochday'] = (d3.index - idate)/pd.Timedelta('1D')

                        d4 = density.orbit_mean(lats=(-61,-29),updown=updown)
                        d4['rrho400'] = 100*(d4['rho400'] - d4['rho400'].mean())/d4['rho400'].mean()
                        d4['epochday'] = (d4.index - idate)/pd.Timedelta('1D')

                        d5 = density.orbit_mean(lats=(-85,-59),updown=updown)
                        d5['rrho400'] = 100*(d5['rho400'] - d5['rho400'].mean())/d5['rho400'].mean()
                        d5['epochday'] = (d5.index - idate)/pd.Timedelta('1D')

                        data[k1] = data[k1].append(d)
                        data1[k1] = data1[k1].append(d1)
                        data2[k1] = data2[k1].append(d2)
                        data3[k1] = data3[k1].append(d3)
                        data4[k1] = data4[k1].append(d4)
                        data5[k1] = data5[k1].append(d5)
        pd.to_pickle((data,data1,data2,data3,data4,data5),'/data/tmp/t6.dat')
        gc.collect()
    def f14():
        """ imf and indices during epoch days of sblist
        Output:
        (imf,index): index=[dataframe,dataframe],in which
        index[0] for away-toward, index[1] for toward-away
        """
        import gc
        lrday = 5
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2011-1-1']
        dataimf = [pd.DataFrame(),pd.DataFrame()]
        dataindex = [pd.DataFrame(),pd.DataFrame()]
        for k1, k2 in enumerate(['away-toward', 'toward-away']):
            sbtmp = sblist[sblist.sbtype==k2]
            for idate in sbtmp.index:
                datelist = idate + pd.TimedeltaIndex(np.arange(10)-lrday, 'd')
                imf = get_imf(datelist)
                sgmindex = get_index(datelist)
                imf['epochday'] = (imf.index - idate)/pd.Timedelta('1D')
                sgmindex['epochday'] = (sgmindex.index - idate)/pd.Timedelta('1D')
                dataimf[k1] = dataimf[k1].append(imf)
                dataindex[k1] = dataindex[k1].append(sgmindex)
        pd.to_pickle((dataimf,dataindex),'/data/tmp/t7.dat')
        gc.collect()
        return
    def f15():
        """ imf and indices variation during epoch days of sblist
        """
        imf, sgmindex = pd.read_pickle('/data/tmp/t7.dat')
        imfgroup = [imf[k].groupby([imf[k].index.month,np.floor(imf[k].epochday*24)])
                for k in [0,1]]
        imfgroup = [imfgroup[k]['Bx','Bye','Bzm'].median() for k in [0,1]]
        indexgroup = [sgmindex[k].groupby(
                [sgmindex[k].index.month, np.floor(sgmindex[k].epochday*24)])
                for k in [0,1]]
        indexgroup = [indexgroup[k]['ap','f107'].median() for k in [0,1]]
        for k in [0,1]:
            imfgroup[k].index.names = ('month', 'epochhour')
            imfgroup[k] = imfgroup[k].reset_index().pivot(index='epochhour', columns='month')
            indexgroup[k].index.names = ('month', 'epochhour')
            indexgroup[k] = indexgroup[k].reset_index().pivot(index='epochhour', columns='month')
        fig,ax = plt.subplots(4,2,sharex=True,sharey=True,figsize=(8,7))
        for k3 in range(2):
            plt.sca(ax[3,k3])
            data = indexgroup[k3]['ap']
            hc1 = plt.contourf(data.columns, data.index/24, data.values,
                               levels=np.linspace(0,20,11))
            for k1,k2 in enumerate(['Bx','Bye','Bzm']):
                plt.sca(ax[k1,k3])
                data = imfgroup[k3][k2]
                hc2 = plt.contourf(data.columns, data.index/24, data.values,
                                   levels=np.linspace(-4,4,11))
        for k1 in range(4):
            for k2 in range(2):
                plt.sca(ax[k1,k2])
                plt.xlim([1,12])
                plt.xticks(np.arange(1,13))
                plt.ylim([-5,5])
                plt.yticks(np.arange(-4,5,2))
        ax[0,0].set_title('away-toward')
        ax[0,1].set_title('toward-away')
        cax1 = plt.axes([0.92,0.4,0.02,0.4])
        plt.colorbar(hc2,cax=cax1,ticks=np.arange(-4,5,2))
        cax2 = plt.axes([0.92,0.1,0.02,0.18])
        plt.colorbar(hc1,cax=cax2,ticks=np.arange(0,21,4))
        plt.show()
        ax[0,1].set_ylabel('Bx',fontsize=14)
        ax[1,1].set_ylabel('By',fontsize=14)
        ax[2,1].set_ylabel('GSM Bz',fontsize=14)
        ax[3,1].set_ylabel('ap',fontsize=14)
        ax[3,0].set_xlabel('Month',fontsize=14)
        ax[3,1].set_xlabel('Month',fontsize=14)
        plt.sca(ax[2,0])
        plt.text(-0.2,1,'Epoch day',fontsize=14,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=plt.gca().transAxes,
                 rotation='vertical')
        return imfgroup
    def f16():
        """ CHAMP and GRACE density variations versus date and epoch days
        """
        rho = list(pd.read_pickle('/data/tmp/t6.dat'))
        fig,ax = plt.subplots(1,2,sharex=False,sharey=True,figsize=(8,4))
        for k in range(2):
            plt.sca(ax[k])
            plt.scatter(rho[0][k].index.dayofyear +
                        rho[0][k].index.hour/24+
                        rho[0][k].index.minute/24/60, rho[0][k].epochday,
                        s=1,linewidth=0,color='k')
            plt.xlim([0,367])
            plt.xticks(np.arange(0,367,60))
            plt.xlabel('Day of year',fontsize=14)
            plt.ylabel('Epoch day',fontsize=14)
        plt.tight_layout()

        fig, ax = plt.subplots(3,2,figsize=(7.5,7.1),sharex=True, sharey=True)
        fig1, ax1 = plt.subplots(3,2,figsize=(7.5,7.1),sharex=True, sharey=True)
        a = [1,1,1,1]
        for k0,k1 in enumerate(['60~90','-30~30','-90~-60']):
            density = rho[k0*2+1]
            for k3,k2 in enumerate(['away-toward','toward-away']):
                dentmp = density[k3]
                if k3==0:
                    dentmp = dentmp[(dentmp.index.month>=2) & (dentmp.index.month<=4)]
                else:
                    dentmp = dentmp[(dentmp.index.month>=8) & (dentmp.index.month<=10)]
                def percentile(n):
                    def percentile_(x):
                        return np.percentile(x,n)
                    percentile_.__name__ = 'percentile_%s' % n
                    return percentile_
                for k4, k5, k6 ,k7, k8, k9 in zip(['gray','red','blue','lightgreen'],
                                                  ['black','darkred','darkblue','darkgreen'],
                                                  [4,3,2,1],
                                                  [1,0.9,0.8,0.7],
                                                  [1,2,3,4],
                                                  ['dawn','noon','dusk','night']):

                    if k8==1:
                        dentmp1 = dentmp[(dentmp.LT>=3) & (dentmp.LT<=9)]
                    elif k8==2:
                        dentmp1 = dentmp[(dentmp.LT>=9) & (dentmp.LT<=15)]
                    elif k8==3:
                        dentmp1 = dentmp[(dentmp.LT>=15) & (dentmp.LT<=21)]
                    else:
                        dentmp1 = dentmp[(dentmp.LT>=21) | (dentmp.LT<=3)]
                    plt.sca(ax[k0,k3])
                    #plt.scatter(dentmp1.epochday,dentmp1.rrho400,
                    #            linewidths=0,color=k4,alpha=k7, zorder=k6)
                    epochhourbin = dentmp1.epochday*24//3*3+1.5
                    dentmp11 = dentmp1.groupby(epochhourbin)['rrho400'].agg(
                            [np.median, percentile(25),percentile(75)])
                    dentmp11.columns = ['median', 'p25', 'p75']
                    a[k8-1], = plt.plot(dentmp11.index/24, dentmp11['median'],zorder=k6+10,color=k5,label=k9)
                    #plt.errorbar(dentmp11.index/24, dentmp11['median'],
                    #             yerr=[dentmp11['median']-dentmp11.p25,dentmp11.p75-dentmp11['median']],
                    #             linewidth=1,zorder=k6+10,color=k5)
                plt.xlim(-5,5)
                plt.xticks(np.arange(-5,6))
                plt.ylim(-40,40)
                plt.yticks(np.arange(-40,41,20))
                plt.grid(which='both')
                if k3 ==0:
                    plt.ylabel(r'$\Delta \rho$ (%)',fontsize=14)
                if k0 ==2:
                    plt.xlabel('Epoch Day',fontsize=14)
                if k0==0 and k3==0:
                    plt.title('Away-Toward')
                if k0==0 and k3==1:
                    plt.title('Toward-Away')

                    #plt.sca(ax1[k0,k3])
                    #plt.scatter(dentmp1.epochday,dentmp1.rho400,
                    #            linewidths=0,color=k4,alpha=k7, zorder=k6)
                    #epochhourbin = dentmp1.epochday*24//3*3+1.5
                    #dentmp12 = dentmp1.groupby(epochhourbin)['rho400'].agg(
                    #        [np.median, percentile(25),percentile(75)])
                    #dentmp12.columns = ['median', 'p25', 'p75']
                    #plt.plot(dentmp12.index/24, dentmp12['median'],zorder=k6+10,color=k5,label=k9)
                    #plt.errorbar(dentmp12.index/24, dentmp12['median'],
                    #             yerr=[dentmp12['median']-dentmp12.p25,dentmp12.p75-dentmp12['median']],
                    #             linewidth=1,zorder=k6+10,color=k5,label=k9)
                    #plt.xlim(-5,5)
                    #plt.xticks(np.arange(-5,6))
                    #plt.ylim(0,6e-12)
                    #plt.yticks(np.arange(0,7,2)*1e-12)
        plt.figure(fig.number)
        plt.legend((a[0],a[1],a[2],a[3]),('Dawn', 'Noon','Dusk','Night'),
                   loc=[-1.3,-0.4],frameon=False,ncol=4)

                #dentmp2 = dentmp[(dentmp.LT>=9) & (dentmp.LT<=15)]
                #plt.scatter(dentmp2.epochday,dentmp2.rrho400,
                #            linewidths=0,color='blue',alpha=0.8, zorder=1)
                #epochhourbin = dentmp2.epochday*24//3*3+1.5
                #dentmp2 = dentmp2.groupby(epochhourbin)['rrho400'].median()
                #plt.plot(dentmp2.index/24, dentmp2,zorder=11,color='darkblue')

                #dentmp3 = dentmp[(dentmp.LT>=15) & (dentmp.LT<=21)]
                #plt.scatter(dentmp3.epochday,dentmp3.rrho400,
                #            linewidths=0,color='red',alpha=0.6,zorder=2)
                #epochhourbin = dentmp3.epochday*24//3*3+1.5
                #dentmp3 = dentmp3.groupby(epochhourbin)['rrho400'].median()
                #plt.plot(dentmp3.index/24, dentmp3,zorder=12,color='darkred')

                #dentmp4 = dentmp[(dentmp.LT>=21) | (dentmp.LT<=3)]
                #plt.scatter(dentmp4.epochday,dentmp4.rrho400,
                #            linewidths=0,color='lightgreen',alpha=1,zorder=3)
                #epochhourbin = dentmp4.epochday*24//3*3+1.5
                #dentmp4 = dentmp4.groupby(epochhourbin)['rrho400'].median()
                #plt.plot(dentmp4.index/24, dentmp4,zorder=14,color='darkgreen')
        return dentmp


        fig,ax = plt.subplots(6,2,sharex=True,sharey=True,figsize=(8,10))
        fig1,ax1 = plt.subplots(3,2,sharex=True,sharey=True,figsize=(8,6))
        for k0,k1 in enumerate(['all', '60~90', '30~60','-30~30','-60~-30','-90~-60']):
            rho[k0] = [rho[k0][k].groupby([rho[k0][k].index.month, np.floor(rho[k0][k].epochday*16)])
                   for k in range(2)]
            rho[k0] = [rho[k0][k]['rrho400'].median() for k in range(2)]
            for k in range(2):
                rho[k0][k].index.names = ['month','epochhour']
                rho[k0][k] = rho[k0][k].reset_index().pivot('epochhour','month')

            for k in range(2):
                plt.sca(ax[k0,k])
                data = rho[k0][k]['rrho400']
                hc = plt.contourf(data.columns, data.index/16, data.values,
                                  levels=np.linspace(-30,30,21))
                plt.xlim([1,12])
                plt.xticks(range(1,13))
                if k is 1:
                    plt.text(1.1,0.5,'Lat: '+k1, horizontalalignment='center',
                            verticalalignment='center',
                            transform=plt.gca().transAxes,
                            rotation='vertical')

                if k0 in [1,3,5]:  # For northern pole, equator, southern pole
                    plt.sca(ax1[np.int((k0-1)/2),k])
                    kk = [1,2,3,4] if k==0 else [7,8,9,10]
                    for k2 in kk:
                        plt.plot(data.index/16, data[k2], color='gray')
                    plt.plot(data.index/16, data[kk].median(axis=1), color='r')
                    plt.xlim(-5,5)
                    plt.ylim(-40,40)
                    plt.xticks(range(-5,6))
                    plt.yticks(np.arange(-40,41,20))
                    if k==0:
                        plt.ylabel(r'$\rho$ (%)')
                    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
                    plt.grid()
        ax[0,0].set_title('away-toward')
        ax[0,1].set_title('toward-away')
        ax[-1,0].set_xlabel('Month',fontsize=14)
        ax[-1,1].set_xlabel('Month',fontsize=14)
        ax1[0,0].set_title('away-toward')
        ax1[0,1].set_title('toward-away')
        ax1[-1,0].set_xlabel('Epoch day',fontsize=14)
        ax1[-1,1].set_xlabel('Epoch day',fontsize=14)
        cax = fig.add_axes([0.3,0.03,0.4,0.01])
        cbar = plt.colorbar(
                hc,cax=cax,
                orientation='horizontal',
                ticks=np.arange(-30,31,10))
        cbar.ax.set_title(r'$\rho$ (%)')
        plt.sca(ax[3,0])
        plt.text(-0.2,1,'Epoch day', horizontalalignment='center',
                verticalalignment='center',
                transform=plt.gca().transAxes,
                rotation='vertical')
        plt.subplots_adjust(bottom=0.12, wspace=0.1)
        plt.show()
        return


    def f17():
        """ MLT/Mlat and LT/Lat distributions of the CHAMP satellites.
        """
        a = get_density_dates(pd.date_range('2005-3-1','2005-3-10'),satellite='champ')
        fig, ax = plt.subplots(2,1,sharex=True)
        plt.sca(ax[0])
        a.add_updown(whichlat='Mlat')
        plt.scatter(a.loc[a.isup,'MLT'], a.loc[a.isup,'Mlat'], linewidths=0,label='up')
        plt.scatter(a.loc[a.isdown,'MLT'], a.loc[a.isdown,'Mlat'], linewidths=0,c='r',label='down')
        plt.xlim(0,24)
        plt.ylim(-95,95)
        plt.xticks(range(0,25,4))
        plt.yticks(range(-90,91,30))
        plt.legend(frameon=False,loc='center left')
        plt.xlabel('MLT')
        plt.ylabel('MLat')

        plt.sca(ax[1])
        a.add_updown()
        plt.scatter(a.loc[a.isup,'LT'], a.loc[a.isup,'lat'], linewidths=0,label='up')
        plt.scatter(a.loc[a.isdown,'LT'], a.loc[a.isdown,'lat'], linewidths=0,c='r',label='down')
        plt.xlim(0,24)
        plt.ylim(-95,95)
        plt.xticks(range(0,25,4))
        plt.yticks(range(-90,91,30))
        plt.legend(frameon=False,loc='center left')
        plt.xlabel('LT')
        plt.ylabel('Lat')
        plt.show()

        fig = plt.figure()
        dt = a.index.dayofyear+a.index.hour/24+a.index.minute/60/24
        h = plt.subplot(polar=True)
        theta = a['MLT']/12*np.pi
        r = 90 + a['Mlat']
        hcup = plt.scatter(theta, r, linewidths=0,c=dt,alpha=0.5)
        set_lt_lat_polar(h,pole='S')
        h.set_rmax(30)
        plt.show()
        return a


    def f18():
        """ Lat/Lon variations
        """
        lday, rday = 5,5
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        rho = [[pd.DataFrame(), pd.DataFrame(), pd.DataFrame()],
               [pd.DataFrame(), pd.DataFrame(), pd.DataFrame()]]  #NH, EQ, SH
        nn = [0,0]
        def ff(x):  # calculate relative value
            x['rrho400'] = 100*(x['rho400']-x['rho400'].mean())/x['rho400'].mean()
            return x
        if False:
            for k0, k in enumerate(['away-toward', 'toward-away']):
                sbtmp = sblist[sblist.sbtype==k]
                for k1 in sbtmp.index:
                    epoch = k1+pd.TimedeltaIndex(np.arange(-lday,rday),'D')
                    for k2 in ['champ', 'grace']:
                        density = get_density_dates(epoch,satellite=k2)
                        if density.empty:
                            print(k2,' has no data around ',k1)
                            continue
                        density = density.add_index()
                        l1= len(density)
                        density = ChampDensity(density[density.Kp<=40])
                        l2= len(density)
                        if l2<=0.8*l1:
                            print('geomagnetic active time around ',k1 )
                            continue
                        nn[k0] = nn[k0]+1
                        print(k,nn)
                        density = density.add_updown()
                        density['epochday'] = (density.index-k1)/pd.Timedelta('1D')
                        # Due to small area, lats>=87 are binned together
                        density.loc[density.lat3==90,'lat3'] = 87
                        density.loc[density.lat3==-90,'lat3'] = -87
                        #  Lat: 60 ~ 90
                        densityn = density[density.lat3>=50]
                        densityn = densityn.groupby(densityn.isup).apply(ff)
                        densityn = densityn[['lat3','long','LT','epochday','rho400','rrho400']]
                        rho[k0][0] = rho[k0][0].append(densityn)
                        #  Lat: -60 ~ 60
                        densitye = density[(density.lat3<=70) & (density.lat3>=-70)]
                        densitye = densitye.groupby(densitye.isup).apply(ff)
                        densitye = densitye[['lat3','long','LT','epochday','rho400','rrho400']]
                        rho[k0][1] = rho[k0][1].append(densitye)
                        #  Lat: -60 ~ -90
                        densitys = density[density.lat3<=-50]
                        densitys = densitys.groupby(densitys.isup).apply(ff)
                        densitys = densitys[['lat3','long','LT','epochday','rho400','rrho400']]
                        rho[k0][2] = rho[k0][2].append(densitys)
            pd.to_pickle(rho,'/data/tmp/t8.dat')
        """
        density = pd.read_pickle('/data/tmp/t8.dat')
        figax = [1,2,3,4,5,6]
        figax[0] = plt.subplots(2*lday,4,figsize=[5,11])
        figax[1] = plt.subplots(2*lday,4,figsize=[7,11])
        figax[2] = plt.subplots(2*lday,4,figsize=[5,11])
        figax[3] = plt.subplots(2*lday,4,figsize=[5,11])
        figax[4] = plt.subplots(2*lday,4,figsize=[7,11])
        figax[5] = plt.subplots(2*lday,4,figsize=[5,11])
        for k0, k in enumerate(['Away-Toward', 'Toward-Away']):
            rho = density[k0]
            for k1,k2 in enumerate(['60~90','-60~60','-60~-90']):
                rho1 = rho[k1]
                if k0==0:
                    rho1 = rho1[(rho1.index.month>=2) &(rho1.index.month<=4)]
                else:
                    rho1 = rho1[(rho1.index.month>=8) &(rho1.index.month<=10)]
                for k3 in np.arange(-lday,rday):   # epoch days
                    for k4,k5 in enumerate(['Dawn','Noon','Dusk','Night']):
                        if k5 is 'Dawn':
                            rho1tmp = rho1[(rho1.LT>3) & (rho1.LT<9)]
                        elif k5 is 'Noon':
                            rho1tmp = rho1[(rho1.LT>9) & (rho1.LT<15)]
                        elif k5 is 'Dusk':
                            rho1tmp = rho1[(rho1.LT>15) & (rho1.LT<21)]
                        else:
                            rho1tmp = rho1[(rho1.LT>21) | (rho1.LT<3)]
                        plt.sca(figax[k0*3+k1][1][k3+lday,k4])
                        rho2 = rho1tmp[np.floor(rho1tmp.epochday)==k3]
                        mlatbin = rho2.lat3
                        mlonbin = np.floor(rho2.long/45)*45+22.5
                        rho2 = rho2.groupby([mlatbin,mlonbin])['rrho400'].median().reset_index().pivot(
                                index='lat3',columns='long',values='rrho400')
                        rho2[rho2.columns[0]+360] = rho2[rho2.columns[0]]

                        if k2 is '60~90':
                            # Area of lat3>87 is small, there should be only one value.
                            rho2.iloc[rho2.index.argmax()]=np.nanmean(rho2.iloc[rho2.index.argmax()])
                            mp = Basemap(projection='npstere',boundinglat=60,lon_0=0,
                                         resolution='c',round=True)
                        elif k2 is '-60~-90':
                            rho2.iloc[rho2.index.argmin()]=np.nanmean(rho2.iloc[rho2.index.argmin()])
                            mp = Basemap(projection='spstere',boundinglat=-60,lon_0=0,
                                         resolution='c',round=True)
                        else:
                            mp = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=60,
                                        llcrnrlon=-180,urcrnrlon=180,resolution='c')
                        mp.drawparallels(np.arange(-80,81,10.),dashes=(1,1),linewidth=0.5)
                        mp.drawmeridians(np.arange(0,360,60.), dashes=(1,1),linewidth=0.5)

                        mlong,mlatg = np.meshgrid(rho2.columns,rho2.index)
                        cf = mp.contourf(mlong,mlatg,rho2.values,
                                         latlon=True,levels=np.linspace(-30,30,21),
                                         extend='both')
                        if k2 is '60~90':
                            mp.scatter(-83.32,82.23,latlon=True,marker='x',c='k')
                        elif k2 is '-60~-90':
                            mp.scatter(126.24,-74.18,latlon=True,marker='x',c='k')
                        if k3==-lday:
                            plt.title(k5)
                        if k4==0 and k1!=1:
                            plt.text(-0.3,0.5,k3,transform=plt.gca().transAxes)
                        elif k4==0 and k1==1:
                            plt.text(-0.6,0.5,k3,transform=plt.gca().transAxes)
                plt.text(0.5,0.95,'{:s}: {:s}'.format(k,k2),
                         horizontalalignment='center',transform=plt.gcf().transFigure)
                plt.subplots_adjust(bottom=0.05,hspace=0.01)
                cax = plt.axes([0.20,0.025,0.6,0.01])
                plt.colorbar(cf, cax=cax, ticks=np.arange(-30,31,10),orientation='horizontal')
        """
        density = pd.read_pickle('/data/tmp/t8.dat')
        figax = [1,2,3,4,5,6]
        figax[0] = plt.subplots(2*lday,4,figsize=[5,11],subplot_kw=dict(polar=True))
        figax[1] = plt.subplots(2*lday,4,figsize=[7,11],sharex=True,sharey=True)
        figax[2] = plt.subplots(2*lday,4,figsize=[5,11],subplot_kw=dict(polar=True))
        figax[3] = plt.subplots(2*lday,4,figsize=[5,11],subplot_kw=dict(polar=True))
        figax[4] = plt.subplots(2*lday,4,figsize=[7,11],sharex=True,sharey=True)
        figax[5] = plt.subplots(2*lday,4,figsize=[5,11],subplot_kw=dict(polar=True))
        for k0, k in enumerate(['Away-Toward', 'Toward-Away']):
            rho = density[k0]
            for k1,k2 in enumerate(['60~90','-60~60','-60~-90']):
                rho1 = rho[k1]
                if k0==0:
                    rho1 = rho1[(rho1.index.month>=2) &(rho1.index.month<=4)]
                else:
                    rho1 = rho1[(rho1.index.month>=8) &(rho1.index.month<=10)]
                for k3 in np.arange(-lday,rday):   # epoch days
                    for k4,k5 in enumerate(['Dawn','Noon','Dusk','Night']):
                        if k5 is 'Dawn':
                            rho1tmp = rho1[(rho1.LT>3) & (rho1.LT<9)]
                        elif k5 is 'Noon':
                            rho1tmp = rho1[(rho1.LT>9) & (rho1.LT<15)]
                        elif k5 is 'Dusk':
                            rho1tmp = rho1[(rho1.LT>15) & (rho1.LT<21)]
                        else:
                            rho1tmp = rho1[(rho1.LT>21) | (rho1.LT<3)]
                        plt.sca(figax[k0*3+k1][1][k3+lday,k4])
                        rho2 = rho1tmp[np.floor(rho1tmp.epochday)==k3]
                        latbin = rho2.lat3
                        utbin = pd.Series(rho2.index.hour//3*3+1.5,index=rho2.index,name='UT')
                        rho2 = rho2.groupby([latbin,utbin])['rrho400'].median().reset_index().pivot(
                                index='lat3',columns='UT',values='rrho400')
                        rho2[rho2.columns[0]+24] = rho2[rho2.columns[0]]
                        rho2[rho2.columns[-2]-24] = rho2[rho2.columns[-2]]
                        rho2 = rho2[np.roll(rho2.columns,1)]

                        if k2 is '60~90':
                            # Area of lat3>87 is small, there should be only one value.
                            rho2.iloc[rho2.index.argmax()]=np.nanmean(rho2.iloc[rho2.index.argmax()])
                            utg,mlatg = np.meshgrid(rho2.columns,rho2.index)
                            thetag, rg = utg/12*np.pi, 90-mlatg
                        elif k2 is '-60~-90':
                            rho2.iloc[rho2.index.argmin()]=np.nanmean(rho2.iloc[rho2.index.argmin()])
                            utg,mlatg = np.meshgrid(rho2.columns,rho2.index)
                            thetag, rg = utg/12*np.pi, 90+mlatg
                        else:
                            utg,mlatg = np.meshgrid(rho2.columns,rho2.index)
                            thetag, rg = utg, mlatg

                        cf = plt.gca().contourf(
                                thetag,rg,rho2.values,
                                levels=np.linspace(-30,30,21),
                                extend='both')
                        if k2 is '-60~60':
                            plt.yticks(np.arange(-60,61,60),fontsize=11)
                            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
                            plt.xlim(0,24)
                            plt.xticks(np.arange(0,25,6))
                        else:
                            plt.gca().set_rmax(30)
                            plt.gca().set_theta_zero_location('S')
                            plt.yticks(np.arange(10,31,20),[''])
                            plt.xticks(np.arange(0,361,60)/180*np.pi,['','','','','',''])
                            if k3 == rday-1 and k4 ==3:
                                plt.gca().set_thetagrids(
                                        np.arange(0,360,60),labels=np.arange(0,24,4),
                                        frac = 1.2,fontsize=10)
                        if k3==-lday:
                            plt.title(k5)
                        if k4==0 and k1!=1:
                            plt.text(-0.3,0.5,k3,transform=plt.gca().transAxes)
                        elif k4==0 and k1==1:
                            plt.text(-0.6,0.5,k3,transform=plt.gca().transAxes)
                plt.text(0.5,0.95,'{:s}: {:s}'.format(k,k2),
                         horizontalalignment='center',transform=plt.gcf().transFigure)
                plt.subplots_adjust(bottom=0.07)
                cax = plt.axes([0.20,0.025,0.6,0.01])
                plt.colorbar(cf, cax=cax, ticks=np.arange(-30,31,10),orientation='horizontal')
        #"""
        return
    def f19():
        """ MLat/MLT variations
        """
        lday, rday = 5,5
        sblist = get_sblist()
        rho = [[pd.DataFrame(), pd.DataFrame()], [pd.DataFrame(),pd.DataFrame()]]
        if False:
            for k0, k in enumerate(['away-toward', 'toward-away']):
                sbtmp = sblist[sblist.sbtype==k]
                for k1 in sbtmp.index:
                    epoch = k1+pd.TimedeltaIndex(np.arange(-lday,rday),'D')
                    density = get_density_dates(epoch)
                    if density.empty:
                        continue
                    density = density.add_index()
                    l1= len(density)
                    density = density[density.Kp<=40]
                    l2= len(density)
                    if l2<=0.7*l1:
                        continue
                    density['epochday'] = (density.index-k1)/pd.Timedelta('1D')
                    densityn = density[density.Mlat>=60]
                    densityn['rrho400'] = 100*(
                            densityn['rho400']-densityn['rho400'].mean())/densityn['rho400'].mean()
                    rho[k0][0] = rho[k0][0].append(densityn)

                    densitys = density[density.Mlat<=-60]
                    densitys['rrho400'] = 100*(
                            densitys['rho400']-densitys['rho400'].mean())/densitys['rho400'].mean()
                    rho[k0][1] = rho[k0][1].append(densitys)
            pd.to_pickle(rho,'/data/tmp/t9.dat')
        density = pd.read_pickle('/data/tmp/t9.dat')
        figax = [plt.subplots(2,lday,figsize=[10,3.5],
                 subplot_kw=dict(projection='polar')) for ddkksfa in range(4)]
        for k0, k in enumerate(['Away-Toward', 'Toward-Away']):
            rho = density[k0]
            for k1,k2 in enumerate(['NH','SH']):
                rho1 = rho[k1]
                if k0==0:
                    rho1 = rho1[(rho1.index.month>=2) &(rho1.index.month<=4)]
                else:
                    rho1 = rho1[(rho1.index.month>=8) &(rho1.index.month<=10)]
                for k3 in np.arange(-lday,rday):   # epoch days
                    plt.sca(figax[k0*2+k1][1][(k3+lday)//lday,(k3+lday)%lday])
                    rho2 = rho1[np.floor(rho1.epochday)==k3]
                    mlatbin = np.floor(rho2.Mlat/3)*3+1.5
                    mltbin = np.floor(rho2.MLT/3)*3+1.5
                    rho2 = rho2.groupby([mlatbin,mltbin])['rrho400'].median().reset_index().pivot(
                            index='Mlat',columns='MLT',values='rrho400')
                    rho2[rho2.columns[0]+24] = rho2[rho2.columns[0]]

                    if k2 is 'NH':
                        rho2.iloc[rho2.index.argmax()]=np.nanmean(rho2.iloc[rho2.index.argmax()])
                        mlt,mlat = np.meshgrid(rho2.columns/12*np.pi,90-rho2.index)
                    else:
                        rho2.iloc[rho2.index.argmin()]=np.nanmean(rho2.iloc[rho2.index.argmin()])
                        mlt,mlat = np.meshgrid(rho2.columns/12*np.pi,90+rho2.index)
                    cf = plt.contourf(mlt,mlat,rho2.values,
                                      levels=np.linspace(-40,40,21),
                                      extend='both')
                    if k2 is 'NH':
                        set_lt_lat_polar(plt.gca(), pole='N',boundinglat=60)
                    else:
                        set_lt_lat_polar(plt.gca(), pole='S',boundinglat=-60)
                cax = plt.axes([0.93,0.2,0.02,0.6])
                plt.colorbar(cf, cax=cax, ticks=np.arange(-30,31,10))
                plt.sca(figax[k0*2+k1][1][0,2])
                plt.title(k+', '+k2)
        return rho2


    def f20():
        """ density variation versus Mlat and MLT during sb epoch days
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        lday, rday = 5, 5
        if False:
            density = [pd.DataFrame(), pd.DataFrame()]  # for at and ta conditions
            for k, k0 in enumerate(['away-toward', 'toward-away']):
                sbtmp = sblist[sblist.sbtype==k0]
                for k1 in sbtmp.index:
                    rho = get_density_dates(
                            k1+pd.TimedeltaIndex(np.arange(-lday, rday),'D'))
                    if rho.empty:
                        print('no data around', k1)
                        continue
                    rho = rho.add_index()
                    l1 = len(rho)
                    rho = rho[rho.Kp<=40]
                    l2 = len(rho)
                    if l2<0.8*l1:
                        print('active geomagnetic condition')
                        continue
                    rho['epochday'] = (rho.index-k1)/pd.Timedelta('1D')
                    mlatbin = rho.Mlat//3*3+3  #  3 degree bin
                    mltbin = rho.MLT//3*3+1.5   #   3 hour bin
                    # Near the poles, the data is sparse
                    # North
                    fp = (mlatbin>=87)
                    mlatbin[fp] = 90
                    mltbin[fp] = 0
                    fp = (mlatbin<87) & (mlatbin>=84)
                    fp1 = (rho.MLT>=6) & (rho.MLT<=18)
                    mltbin[(fp) & (fp1)] = 12
                    mltbin[(fp) & (~fp1)] = 0
                    fp = (mlatbin<84) & (mlatbin>=81)
                    mltbin[fp] = rho.loc[fp,'MLT']//6*6+3
                    # South
                    fp = (mlatbin<=-87)
                    mlatbin[fp] = -90
                    mltbin[fp] = 0
                    fp = (mlatbin>-87) & (mlatbin<=-84)
                    fp1 = (rho.MLT>=6) & (rho.MLT<=18)
                    mltbin[(fp) & (fp1)] = 12
                    mltbin[(fp) & (~fp1)] = 0
                    fp = (mlatbin>-84) & (mlatbin<=-81)
                    mltbin[fp] = rho.loc[fp,'MLT']//6*6+3

                    rho = rho.groupby([mlatbin,mltbin])
                    rho = rho.filter(lambda x: len(np.unique(np.floor(x['epochday'])))>=lday+rday)
                    if rho.empty:
                        print('There is data gap around', k1)
                        continue
                    else:
                        print('OK')
                    rho = rho.groupby([mlatbin,mltbin])
                    def ff(x):
                        x['rrho400'] = 100*(x['rho400']-x['rho400'].mean())/x['rho400'].mean()
                        return x
                    rho = rho.apply(ff)
                    rho = rho[['Mlat','MLT', 'epochday', 'rrho400']]
                    density[k] = density[k].append(rho)
            pd.to_pickle(density, '/data/tmp/t10.dat')
        density = pd.read_pickle('/data/tmp/t10.dat')
        fig, ax = plt.subplots(10,4,figsize=[6,10])
        for k, k0 in enumerate(['away-toward', 'toward-away']):
            rho = density[k]
            if k==0:
                fp = (rho.index.month>=2) & (rho.index.month<=4)
                rho = rho[fp]
            else:
                fp = (rho.index.month>=8) & (rho.index.month<=10)
                rho = rho[fp]
            mlatbin = rho.Mlat//3*3+1.5  #  3 degree bin
            mltbin = rho.MLT//3*3+1.5   #  3 hour bin
            # Near the poles, the data is sparse
            # North
            fp = (mlatbin>=87)
            mlatbin[fp] = 90
            mltbin[fp] = 0
            fp = (mlatbin<87) & (mlatbin>=84)
            fp1 = (rho.MLT>=6) & (rho.MLT<=18)
            mltbin[(fp) & (fp1)] = 12
            mltbin[(fp) & (~fp1)] = 0
            fp = (mlatbin<84) & (mlatbin>=81)
            mltbin[fp] = rho.loc[fp,'MLT']//6*6+3
            # South
            fp = (mlatbin<=-87)
            mlatbin[fp] = -90
            mltbin[fp] = 0
            fp = (mlatbin>-87) & (mlatbin<=-84)
            fp1 = (rho.MLT>=6) & (rho.MLT<=18)
            mltbin[(fp) & (fp1)] = 12
            mltbin[(fp) & (~fp1)] = 0
            fp = (mlatbin>-84) & (mlatbin<=-81)
            mltbin[fp] = rho.loc[fp,'MLT']//6*6+3

            epochbin = rho.epochday//1  # 1 day bin
            rho = rho.groupby([epochbin,mlatbin,mltbin])['rrho400'].mean().reset_index([1,2])
            for  k1 in range(-lday, rday):
                for k2,k3 in enumerate(['NH','SH']):
                    plt.sca(ax[k1+lday,k*2+k2])
                    rhotmp = rho[rho.index==k1]
                    if k2==0:
                        rhotmp = rhotmp[rhotmp.Mlat>=0]# Northern hemisphere
                        r = 90-rhotmp.Mlat
                    else:
                        rhotmp = rhotmp[rhotmp.Mlat<=-0]# Northern hemisphere
                        r = 90+rhotmp.Mlat
                    theta = rhotmp.MLT/12*np.pi+np.pi/2
                    x = r*np.cos(theta)
                    y = r*np.sin(theta)
                    z = np.array(rhotmp.rrho400)
                    xg, yg = np.meshgrid(np.arange(-90,91,3),np.arange(-90,91,3))
                    zg = griddata((x,y),z,(xg,yg))
                    plt.contourf(xg,yg,zg,levels=np.arange(-40,41,3))
                    plt.axis('equal')
                    plt.gca().set_xticks([])
                    plt.gca().set_yticks([])
                    plt.gca().spines['right'].set_color('none')
                    plt.gca().spines['top'].set_color('none')
                    plt.gca().spines['bottom'].set_color('none')
                    plt.gca().spines['left'].set_color('none')
                    plt.xlim(-90,90)
                    plt.ylim(-90,90)
        return rho

    def f21():
        """ density variation versus Mlat and Mlon during sb epoch days
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        lday, rday = 5, 5
        if False:
            density = [pd.DataFrame(), pd.DataFrame()]  # for at and ta conditions
            for k, k0 in enumerate(['away-toward', 'toward-away']):
                sbtmp = sblist[sblist.sbtype==k0]
                for k1 in sbtmp.index:
                    rho = get_density_dates(
                            k1+pd.TimedeltaIndex(np.arange(-lday, rday),'D'))
                    if rho.empty:
                        print('no data around', k1)
                        continue
                    rho = rho.add_index()
                    l1 = len(rho)
                    rho = rho[rho.Kp<=40]
                    l2 = len(rho)
                    if l2<0.8*l1:
                        print('active geomagnetic condition')
                        continue
                    rho['epochday'] = (rho.index-k1)/pd.Timedelta('1D')
                    latbin = rho.lat//3*3+3  #  3 degree bin
                    lonbin = rho.long//45*45+22.5   #   3 hour bin
                    # Near the poles, the data is sparse
                    # North
                    fp = (latbin>=87)
                    latbin[fp] = 90
                    lonbin[fp] = 0
                    fp = (latbin<87) & (latbin>=84)
                    fp1 = (rho.long>=0)
                    lonbin[(fp) & (fp1)] = 90
                    lonbin[(fp) & (~fp1)] = -90
                    fp = (latbin<84) & (latbin>=81)
                    lonbin[fp] = rho.loc[fp,'long']//90*90+45
                    # South
                    fp = (latbin<=-87)
                    latbin[fp] = -90
                    lonbin[fp] = 0
                    fp = (latbin>-87) & (latbin<=-84)
                    fp1 = (rho.long>=0)
                    lonbin[(fp) & (fp1)] = 90
                    lonbin[(fp) & (~fp1)] = -90
                    fp = (latbin>-84) & (latbin<=-81)
                    lonbin[fp] = rho.loc[fp,'long']//90*90+45

                    rho = rho.groupby([latbin,lonbin])
                    rho = rho.filter(lambda x: len(np.unique(np.floor(x['epochday'])))>=lday+rday)
                    if rho.empty:
                        print('There is data gap around', k1)
                        continue
                    else:
                        print('OK')
                    rho = rho.groupby([latbin,lonbin])
                    def ff(x):
                        x['rrho400'] = 100*(x['rho400']-x['rho400'].mean())/x['rho400'].mean()
                        return x
                    rho = rho.apply(ff)
                    rho = rho[['lat','long', 'epochday', 'rrho400']]
                    density[k] = density[k].append(rho)
            pd.to_pickle(density, '/data/tmp/t11.dat')
        density = pd.read_pickle('/data/tmp/t11.dat')
        fig, ax = plt.subplots(10,4,figsize=[6,10])
        for k, k0 in enumerate(['away-toward', 'toward-away']):
            rho = density[k]
            if k==0:
                fp = (rho.index.month>=2) & (rho.index.month<=4)
                rho = rho[fp]
            else:
                fp = (rho.index.month>=8) & (rho.index.month<=10)
                rho = rho[fp]
            latbin = rho.lat//3*3+1.5  #  3 degree bin
            lonbin = rho.long//45*45+22.5   #  3 hour bin
            # Near the poles, the data is sparse
            # North
            fp = (latbin>=87)
            latbin[fp] = 90
            lonbin[fp] = 0
            fp = (latbin<87) & (latbin>=84)
            fp1 = (rho.long>=0)
            lonbin[(fp) & (fp1)] = 90
            lonbin[(fp) & (~fp1)] = -90
            fp = (latbin<84) & (latbin>=81)
            lonbin[fp] = rho.loc[fp,'long']//90*90+45
            # South
            fp = (latbin<=-87)
            latbin[fp] = -90
            lonbin[fp] = 0
            fp = (latbin>-87) & (latbin<=-84)
            fp1 = (rho.long>=0)
            lonbin[(fp) & (fp1)] = 90
            lonbin[(fp) & (~fp1)] = -90
            fp = (latbin>-84) & (latbin<=-81)
            lonbin[fp] = rho.loc[fp,'long']//90*90+45

            epochbin = rho.epochday//1  # 1 day bin
            rho = rho.groupby([epochbin,latbin,lonbin])['rrho400'].mean().reset_index([1,2])
            for  k1 in range(-lday, rday):
                for k2,k3 in enumerate(['NH','SH']):
                    plt.sca(ax[k1+lday,k*2+k2])
                    rhotmp = rho[rho.index==k1]
                    if k2==0:
                        rhotmp = rhotmp[rhotmp.lat>=0]# Northern hemisphere
                        r = 90-rhotmp.lat
                    else:
                        rhotmp = rhotmp[rhotmp.lat<=-0]# Northern hemisphere
                        r = 90+rhotmp.lat
                    theta = rhotmp.long/180*np.pi
                    x = r*np.cos(theta)
                    y = r*np.sin(theta)
                    z = np.array(rhotmp.rrho400)
                    xg, yg = np.meshgrid(np.arange(-90,91,3),np.arange(-90,91,3))
                    zg = griddata((x,y),z,(xg,yg))
                    plt.contourf(xg,yg,zg,levels=np.arange(-30,31,3))
                    plt.axis('equal')
                    plt.gca().set_xticks([])
                    plt.gca().set_yticks([])
                    plt.gca().spines['right'].set_color('none')
                    plt.gca().spines['top'].set_color('none')
                    plt.gca().spines['bottom'].set_color('none')
                    plt.gca().spines['left'].set_color('none')
                    plt.xlim(-90,90)
                    plt.ylim(-90,90)
        return rho


    def f22():
        """ density variation for different Mlts during sb epoch days
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        lday, rday = 5, 5
        if False:
            density = [pd.DataFrame(), pd.DataFrame()]  # for at and ta conditions
            for k, k0 in enumerate(['away-toward', 'toward-away']):
                sbtmp = sblist[sblist.sbtype==k0]
                for k1 in sbtmp.index:
                    rho = get_density_dates(
                            k1+pd.TimedeltaIndex(np.arange(-lday, rday),'D'))
                    if rho.empty:
                        print('no data around', k1)
                        continue
                    rho = rho.add_index()
                    l1 = len(rho)
                    rho = rho[rho.Kp<=40]
                    l2 = len(rho)
                    if l2<0.8*l1:
                        print('active geomagnetic condition')
                        continue
                    rho['epochday'] = (rho.index-k1)/pd.Timedelta('1D')
                    mlatbin = rho.Mlat//30*30+15  #  3 degree bin
                    tmp = rho.loc[(rho.MLT>3) & (rho.MLT<9),'rho400']
                    rho.loc[(rho.MLT>3) & (rho.MLT<9),'rrho400']=100*(tmp-tmp.mean())/tmp.mean()
                    tmp = rho.loc[(rho.MLT>9) & (rho.MLT<15),'rho400']
                    rho.loc[(rho.MLT>9) & (rho.MLT<15),'rrho400']=100*(tmp-tmp.mean())/tmp.mean()
                    tmp = rho.loc[(rho.MLT>15) & (rho.MLT<21),'rho400']
                    rho.loc[(rho.MLT>15) & (rho.MLT<21),'rrho400']=100*(tmp-tmp.mean())/tmp.mean()
                    tmp = rho.loc[(rho.MLT>21) | (rho.MLT<3),'rho400']
                    rho.loc[(rho.MLT>21) | (rho.MLT<3),'rrho400']=100*(tmp-tmp.mean())/tmp.mean()
                    rho = rho[['Mlat','MLT', 'epochday', 'rho400','rrho400']]
                    density[k] = density[k].append(rho)
            pd.to_pickle(density, '/data/tmp/t12.dat')
        density = pd.read_pickle('/data/tmp/t12.dat')
        fig, ax = plt.subplots(4,2,sharex=True,sharey=True)
        for k, k0 in enumerate(['away-toward', 'toward-away']):
            rho = density[k]
            if k==0:
                fp = (rho.index.month>=2) & (rho.index.month<=4)
                rho = rho[fp]
            else:
                fp = (rho.index.month>=8) & (rho.index.month<=10)
                rho = rho[fp]
            rho = rho[rho.Mlat>=60]
            for k1,k2 in enumerate(['dawn','noon','dusk','night']):
                plt.sca(ax[k1,k])
                if k1==0:
                    rhotmp=rho[(rho.MLT>3) & (rho.MLT<9)]
                elif k1==1:
                    rhotmp=rho[(rho.MLT>9) & (rho.MLT<15)]
                elif k1==2:
                    rhotmp=rho[(rho.MLT>15) & (rho.MLT<21)]
                else:
                    rhotmp=rho[(rho.MLT>21) | (rho.MLT<3)]
                epochbin = rhotmp.epochday*24//3*3
                rhoz = rhotmp.groupby(epochbin)['rrho400'].mean()
                plt.plot(rhoz.index/24,rhoz)
                plt.ylim(-40,60)
                plt.xlim(-5,5)
                plt.xticks(np.arange(-5,6))
                plt.yticks(np.arange(-40,61,20))
                plt.grid()
    def f23():
        """ For two special cases: 02-10-11~02-10-13 and 08-10-14~08-10-16
        """
        dates = [pd.date_range('2002-10-11','2002-10-13'),
                 pd.date_range('2008-10-14','2008-10-16')]
        fig1,ax1 = plt.subplots(4,2,sharex=True,sharey=True)
        fig2,ax2 = plt.subplots(4,2,sharex=True,sharey=True)
        for k00,k in enumerate(dates):
            for k11,k1 in enumerate(['champ','grace']):
                rho = get_density_dates(k,k1)
                rho = rho.add_updown()
                idx = get_index(k)
                imf = get_imf(k)
                for k22,k2 in enumerate(['up','down']):
                    if k2=='up':
                        rho1 = rho[rho.isup]
                    if k2=='down':
                        rho1 = rho[rho.isdown]
                    for k3 in np.arange(-72,-60,3):
                        if k00==0:
                            plt.sca(ax1[k22*2+1,k11])
                        if k00==1:
                            plt.sca(ax2[k22*2+1,k11])
                        rho2 = rho1[rho1.lat3==k3]
                        if rho2.empty:
                            continue
                        plt.plot(rho2.index, rho2.rho400/1e-12,'k')
                    for k3 in np.arange(60,72,3):
                        if k00==0:
                            plt.sca(ax1[k22*2+0,k11])
                        if k00==1:
                            plt.sca(ax2[k22*2+0,k11])
                        rho2 = rho1[rho1.lat3==k3]
                        if rho2.empty:
                            continue
                        plt.plot(rho2.index, rho2.rho400/1e-12,'k')
    def f24():
        """Near the geodetic poles, the longitude and LT are not important;
        so it is a good location for the research of UT variation
        """
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2010-12-31']
        for k in sblist.index
#--------------------------#
    plt.close('all')
    a = f16()
    plt.show()
    import gc
    gc.collect()
