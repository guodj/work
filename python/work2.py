#-*- coding: utf-8 -*-

"""For the second work
"""

__author__ = 'Dongjie Guo'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import griddata
import pdb   # set breakpoint
import matplotlib.dates as mdates
from scipy.signal import argrelextrema


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
    sblist.ix[(doy>35)  & (doy<125),'season'] = 'me'
    sblist.ix[(doy>221) & (doy<311),'season'] = 'se'
    sblist.ix[(doy>128) & (doy<218),'season'] = 'js'
    sblist.ix[(doy>311) | (doy<36),'season'] = 'ds'
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


def get_cirlist1():
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

    def add_updown(self, latcolumnname='lat'):
        """ add 'isup' and 'isdown' columns to self
        Note that the function is appropriate for continuous data

        Args:
            latcolumnname: data column name used to seperate up and down
                orbits
        Returns:
            self added with columns 'isup' and 'isdown'

        Note:
            The results may be wrong near poles. But they are not excluded.
            Because if doing so, the function should not be run 2 times.
            So please exclude them at the very necessary location.
            Some bugs may exist at the data gap.
        """
        if not self.empty:
            lat = self[latcolumnname]
            dlat = lat.diff()
            dlat.iloc[0] = dlat.iloc[1]
            self['isup'] = (dlat > 0)
            self['isdown'] = (dlat < 0)
            self = self[(self.isup) | (self.isdown)]
            #self = self[np.abs(self.lat3)<self.lat3.max()]
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
            lats: two elements tuple, selected latitudes: lats[0]<=lat<=lats[1]
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


    def relative_density_to_previous_orbit(self, whichcolumn='rho400'):
        """Get the time interval (dhour_po) between two continuous orbits and
        density relative to in the previous orbit (rrho_po).
        (rho-rho_po)/(rho_po*dhour), unit: /h

        Input:
            whichcolumn: 'rho', 'rho400', 'rho410' ...

        Output:
            self with columns 'dhour_po' and 'rrho_po' added
        """
        self = self.add_updown()
        groupdensity = self.groupby(['isup','lat3'])
        def rtl(x):
            x['dhour_po'] = (np.insert(np.diff(x.index),0,'nat')/
                    np.timedelta64(1,'h'))
            x['rrho_po'] = (x[whichcolumn].diff()/
                            np.roll(x[whichcolumn],1)/
                            x['dhour_po'])
            return x
        self = groupdensity.apply(rtl)
        fp2 = (self.dhour_po>1.7)
        fp3 = (self.dhour_po<1.3)
        self.loc[(fp2) | (fp3),'rrho_po'] = np.nan

        #  relative changes for the maximum and minimum latitudes are unreliable
        #  due to algorithm of add_updown
        fp1 = self.lat3<self.lat3.max()
        fp2 = self.lat3>self.lat3.min()
        self = self[(fp1) & (fp2)]
        return ChampDensity(self)


    def difference_density_to_previous_orbit(self, whichcolumn='rho400'):
        """Get the time interval (dhour_po) between two continuous orbits and
        the difference density between the current and  previous orbit, drho_po:
        (rho-rho_po)/dhour, unit: kg/m^3/h

        Input:
            whichcolumn: 'rho', 'rho400', 'rho410' ...

        Output:
            self with columns 'dhour_po' and 'drho_po' added
        """
        self = self.add_updown()
        groupdensity = self.groupby(['isup','lat3'])
        def rtl(x):
            x['dhour_po'] = (np.insert(np.diff(x.index),0,'nat')/
                    np.timedelta64(1,'h'))
            x['drho_po'] = (x[whichcolumn].diff()/x['dhour_po'])
            return x
        self = groupdensity.apply(rtl)
        fp2 = (self.dhour_po>1.7)
        fp3 = (self.dhour_po<1.3)
        self.loc[(fp2) | (fp3),'drho_po'] = np.nan

        #  relative changes for the maximum and minimum latitudes are unreliable
        #  due to algorithm of add_updown
        fp1 = self.lat3<self.lat3.max()
        fp2 = self.lat3>self.lat3.min()
        self = self[(fp1) & (fp2)]
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


def set_lt_lat_polar(h, pole='N'):
    """Change some default sets of a LT-latitude polar coordinates

    Input
    h: axis handle
    pole: 'N' or 'S', Northern or Southern Hemisphere
    """
    h.set_theta_zero_location('S')
    h.set_thetagrids(angles=[0,90,180,270], labels=[0,6,12,18], fontsize=14)
    h.set_rmax(90)
    if pole =='N':
        h.set_rgrids(radii=[30,60,90],
                     labels=['$60^\circ$', '$30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    if pole =='S':
        h.set_rgrids(radii=[30,60,90],
                     labels=['$-60^\circ$', '$-30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)


def set_lon_lat(h, pole='N'):
    """change some default sets of a longitude-latitude polar coordinates

    Input
    h: axis handle
    pole: 'N' or 'S', Northern or Southern Hemisphere
    """
    h.set_theta_zero_location('S')
    ax.set_thetagrids(
            (-90,0,90,180),
            labels=('$-90^\circ$','$0^\circ$', '$90^\circ$','$\pm180^\circ$'),
            fontsize=14)
    h.set_rmax(90)
    if pole =='N':
        h.set_rgrids(radii=[30,60,90],
                     labels=['$60^\circ$', '$30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)
    if pole =='S':
        h.set_rgrids(radii=[30,60,90],
                     labels=['$-60^\circ$', '$-30^\circ$', '$0^\circ$'],
                     angle=135, fontsize=14)


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
        rho = rho.drop_duplicates()
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


if __name__=='__main__':

    def f1(): # test relative_density_to_previous_orbit and contourf_date_lat
        density = get_density_dates(pd.date_range('2003-10-27','2003-11-3'))
        density = density.relative_density_to_previous_orbit()
        f = plt.figure()
        plt.subplot(211)
        h = density.contourf_date_lat(
                plt.gca(),whichcolumn='rrho_po', updown='down',
                cmap='bwr', levels=np.linspace(-0.5,0.5,11),
                extend='both')
        plt.colorbar(h)
        plt.grid(which='both')

        plt.subplot(212)
        h = density.contourf_date_lat(
                plt.gca(),whichcolumn='rho400',updown='down')
        plt.colorbar(h)
        plt.tight_layout()
        plt.show()

    def f2():    # relative density to the 0 epoch time.
        satellite = 'champ'
        fpath = ('/data/tmp/t0.dat')
        if False:      # data preparation
            density = get_density_dates(pd.date_range('2001-1-1','2010-12-31'),
                                        satellite=satellite)
            density = density.relative_density_to_previous_orbit(
                    whichcolumn='rho400')   # updown columns are added here
            density = density.add_epoch_sb(exclude_cir=True)    # add epochat and epochta

            density = density.add_imf()

            #density = density.add_index()   # add ap column
            #density = ChampDensity(density[density.ap<=20]) # geomagnetic active time is excluded
            density = density[['lat3','epochta', 'epochat','rrho_po',
                               'Bx', 'Bym', 'Bzm']]
            density.to_pickle(fpath)
            import gc
            gc.collect()
        density = pd.read_pickle(fpath)
        density = density[(density.lat3>=-84) & (density.lat3<=84)]
        # bin data
        dt = 3/24        # the bin window in the epoch day
        lrday = 5
        density['epochta'] = pd.cut(
                density['epochta'],
                bins=np.arange(-lrday,lrday+1e-6,dt),
                labels=np.arange(-lrday,lrday+1e-6,dt)[:-1] + dt/2,
                include_lowest=True).astype('float')
        density['epochat'] = pd.cut(
                density['epochat'],
                bins=np.arange(-lrday,lrday+1e-6,dt),
                labels=np.arange(-lrday,lrday+1e-6,dt)[:-1] + dt/2,
                include_lowest=True).astype('float')
        density['lat9'] = density.lat3//9*9+3
        density =ChampDensity(density)

        f, ax = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(7,8.7)) # density map
        f1, ax1 = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(7,8.7)) # density
        for k1, season in enumerate(['me', 'se', 'js', 'ds']):
            tmp = density.subset_season(season)  # for september equinox
            for k2, epochcolumn in enumerate(['epochat', 'epochta']):
                tmp_sum = (
                        tmp.groupby([epochcolumn,'lat3'])['rrho_po'].
                        size().
                        reset_index().
                        pivot(index='lat3', columns=epochcolumn,values=0))
                tmp_median = (
                        tmp.groupby([epochcolumn,'lat3'])['rrho_po'].
                        median().
                        reset_index().
                        pivot(index='lat3', columns=epochcolumn,
                              values='rrho_po'))

                tmp_median = tmp_median*dt*24+1
                tmp_ts = tmp_median.loc[:, tmp_median.columns<0].prod(axis=1)
                tmp_cum = tmp_median.cumprod(axis=1)
                tmp1 = (tmp_cum.divide(tmp_ts, axis=0)-1)*100
                #tmp1 = pd.rolling_mean(tmp1, window=6,center=True,min_periods=1)

                plt.sca(ax[k1,k2])
                laty = tmp1.index.values
                epochdayx = tmp1.columns.values
                [epochdayx, laty] = np.meshgrid(epochdayx, laty)

                hp = plt.pcolormesh(epochdayx, laty, tmp1.values,
                                    cmap='jet',vmin=-10, vmax=30)
                plt.xlim([-lrday,lrday])
                plt.ylim([-90,90])
                plt.xticks(np.arange(-lrday,lrday+1))
                plt.yticks(np.arange(-90, 91, 45))

                #  UT variation for high latitudes and equator
                plt.sca(ax1[k1,k2])
                #tmp2 = tmp_median*dt*24+1
                Nhigh = tmp1[(tmp1.index>=60) & (tmp1.index<=90)]
                Nhigh = Nhigh.mean(axis=0)
                equator = tmp1[(tmp1.index>=-30) & (tmp1.index<=30)]
                equator = equator.mean(axis=0)
                Shigh = tmp1[(tmp1.index>=-90) & (tmp1.index<=-60)]
                Shigh = Shigh.mean(axis=0)
                plt.plot(Nhigh, color='r', linewidth=2,
                         label=r'$60^\circ$-$90^\circ$',zorder=2, alpha=0.7)
                plt.plot(equator, color='k', linewidth=2,
                         label=r'$-30^\circ$-$30^\circ$',zorder=0 )
                plt.plot(Shigh, color='b', linewidth=2,
                         label=r'$-90^\circ$-$60^\circ$', zorder=1)
                plt.xlim([-lrday,lrday])
                plt.ylim([-35,50])
                plt.xticks(np.arange(-lrday,lrday+1))
                plt.yticks(np.arange(-30, 46, 15))
                plt.grid(axis='x')
                plt.axhline(color='k', linestyle='--')
        ax[0,0].set_title('Away-Toward')
        ax[0,1].set_title('Toward-Away')
        ax[3,0].set_xlabel('Epoch day', fontsize=14)
        ax[3,1].set_xlabel('Epoch day', fontsize=14)
        for k in np.arange(4):
            ax[k,0].set_ylabel('Latitude',fontsize=14)
        f.subplots_adjust(wspace=0.04,bottom=0.1)
        cax = f.add_axes([0.93,0.2,0.01,0.6])
        f.colorbar(hp, cax=cax)
        #---------------------------------------------#
        ax1[0,0].set_title('Away-Toward')
        ax1[0,1].set_title('Toward-Away')
        ax1[3,0].set_xlabel('Epoch day', fontsize=14)
        ax1[3,1].set_xlabel('Epoch day', fontsize=14)
        for k in np.arange(4):
            ax1[k,0].set_ylabel(r'$\Delta\rho$ (%)',fontsize=14)
        ax1[0,0].legend(ncol=3, fontsize=12, frameon=False,
                        loc=(0.25,1.23))
        f1.subplots_adjust(wspace=0.04,top=0.9)
        plt.show()

    def f3():    # daily average density variation between 2001-2010
        fpath = ('/data/tmp/t1.dat')
        if False:
            density1 = get_density_dates(pd.date_range('2001-1-1','2010-12-31'),
                                        satellite='champ')
            density2 = get_density_dates(pd.date_range('2001-1-1','2010-12-31'),
                                        satellite='grace')
            density1 = density1.resample('1D',how='mean')
            density1 = density1[['rho','rho400','rho410']]
            density2 = density2.resample('1D',how='mean')
            density2 = density2[['rho','rho400','rho410']]
            fpath = ('/home/gdj/work/density_sector_champ_2nd/output/data/tmp'
                     '/tmp.dat')
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
        index = index.resample('1D', how='mean')
        import gc
        gc.collect()
        ax2 = ax.twinx()
        ax2.plot(index['f107'],'o',color='grey',markersize=3,linewidth=1,
                 markeredgewidth=0,zorder=3,alpha=0.8,label='F10.7')
        ax2.set_ylim([-100,300])
        ax2.set_yticklabels([])
        ax2.legend(loc=(0.7,0.74),frameon=False)

        plt.show()

    def f4():    # CHAMP or GRACE LT distributions
        satellite = 'grace'
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2006-12-31']
        cirlist = get_cirlist()
        sblist = sblist[cirlist.min().date():cirlist.max().date()]
        sblist = sblist[~sblist.index.isin(pd.to_datetime(cirlist.date))]
        groupedsblist = sblist.groupby([sblist.season, sblist.sbtype])
        def f(x):
            y = [get_density_dates([k],satellite=satellite).
                 LT_median() for k in x.index]
            return np.array(y)
        lt = groupedsblist.apply(f)

        fig, ax = plt.subplots(4,2,sharex=True,sharey=True,figsize=(6.5,6.9))
        for k1, season in enumerate(['me', 'se', 'js', 'ds']):
            for k2, sbtype in enumerate(['away-toward', 'toward-away']):

                plt.sca(ax[k1,k2])
                tmp = np.concatenate((lt[season,sbtype][:,0],
                                      lt[season,sbtype][:,1]))
                plt.hist(tmp,np.arange(0,25,2),color='gray')
                plt.xlim(0,24)
                plt.xticks(np.arange(0,25,4))
                plt.ylim(0,6)
                plt.yticks(np.arange(0,11,2))
        ax[-1,0].set_xlabel('Local Time', fontsize=14)
        ax[-1,-1].set_xlabel('Local Time', fontsize=14)
        plt.subplots_adjust(left=0.1,right=0.95,wspace=0.1,hspace=0.1)
        plt.figtext(0.03,0.5,'{:s} LT distribution'.format(satellite),
                    verticalalignment='center',
                    rotation='vertical',fontsize=14)
        plt.show()


    #  density response to cir, test relative_density_to_previous_orbit
    #  and difference_density_to_previous_orbit
    def f5():
        satellite = 'champ'
        lepoch = 2
        repoch = 7
        deltat = 1.5/24   #  epoch time bin window

        fpath = ('/data/tmp/t2.dat')
        if False:  # data preparation
            density = get_density_dates(pd.date_range('2008-1-1', '2008-12-31'),
                                        satellite='champ')
            density = density.relative_density_to_previous_orbit('rho400')
            density = density.difference_density_to_previous_orbit('rho400')
            cirlist = get_cirlist1()
            density = density.add_epochtime(cirlist,lday=lepoch,rday=repoch)
            density = ChampDensity(density[~np.isnan(density.epochtime)])  #  exclude nan
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
        plt.xlim(-lepoch,repoch)
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
        plt.xlim(-2,7)
        plt.xlabel('Epoch days', fontsize=14)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,30))
        plt.ylabel('Lat (deg)',fontsize=14)
        plt.colorbar(label=r'$\Delta\rho\ (10^{-13} kg/m^{-3})$')

        plt.tight_layout()
        plt.show()


    # another method for calculating density response to CIR
    # the result shows that this method may be better
    def f6():
        ldays = 3
        rdays = 7
        fpath = '/data/tmp/t3.dat'
        if False:  # data preparation
            cirlist = get_cirlist1()
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
        cirlist = get_cirlist1()
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
        """ find SBs that is not accompanied by CIR
        A SB is selected when no CIR occurs during SB epoch days.
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
        if True:  # Show the result.
            a = sblist[sblist.sbtype=='away-toward'].groupby('season').apply(len)
            plt.bar(np.arange(4)+0.1, a,linewidth=0,color='k')
            plt.xlim(0,4)
            plt.xticks(np.arange(0.5,4.5),a.index)
        return sblist


    def f9():
        """ Superposed epoch analysis
        density variations during 10 sector polarity reversal centered days
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
        """Case study of the density response to SBs without CIR
        """
        sblist = f8()
        print(len(sblist))
        plt.close()
        for k in sblist.index:
            adates = k+pd.TimedeltaIndex(np.arange(-3,3),'D')
            sw = get_sw(adates)
            imf = get_imf(adates)
            sgindex = get_index(adates)
            rho = get_density_dates(adates)
            rhonup = rho.orbit_mean(lats=[60,90],updown='up')
            LTnup = rhonup.LT.median()
            rhondown = rho.orbit_mean(lats=[60,90],updown='down')
            LTndown = rhondown.LT.median()

            rhosup = rho.orbit_mean(lats=[-90,-60],updown='up')
            LTsup = rhosup.LT.median()
            rhosdown = rho.orbit_mean(lats=[-90,-60],updown='down')
            LTsdown = rhosdown.LT.median()
            fig,ax = plt.subplots(7,1,sharex=True,figsize=(7,10))
            ax[0].plot(sw.index, sw.speed,'darkgreen',linewidth=2)
            ax[0].set_ylabel(r'$\mathsf{V_{SW}}$',fontsize=14)

            ax[1].plot(imf.index, imf.Bx,'darkblue',linewidth=2)
            ax[1].plot(imf.index, imf.Bye,'darkred',linewidth=2)
            ax[1].axhline(0,color='gray',linestyle='--',zorder=0)
            ax[1].set_ylabel(r'$\mathsf{B_x,\ B_y}$',fontsize=14)

            ax[2].plot(imf.index, imf.Bzm,'darkcyan',linewidth=2)
            ax[2].axhline(0,color='gray',linestyle='--',zorder=0)
            ax[2].set_ylabel(r'$\mathsf{B_z}$',fontsize=14)

            ax[3].plot(sgindex.index, sgindex.f107,'darkorange',linewidth=2)
            ax[3].set_ylabel(r'$\mathsf{F10.7}$',fontsize=14)

            ax[4].plot(sgindex.index, sgindex.ap,'purple',linewidth=2)
            ax[4].set_ylabel(r'$\mathsf{ap}$',fontsize=14)

            ax[5].plot(rhonup.index, rhonup.rho400/1e-12,'gray',
                       linewidth=2,label='LT: {:.1f}'.format(LTnup),zorder=1)
            ax[5].plot(rhondown.index, rhondown.rho400/1e-12,'k',
                       linewidth=2,label='LT: {:.1f}'.format(LTndown),zorder=0)
            ax[5].legend(frameon=False,loc=0)
            ax[5].set_ylabel(r'$\mathsf{\rho,\ 400 km}$',fontsize=14)
            ax[6].plot(rhosup.index, rhosup.rho400/1e-12,'gray',
                       linewidth=2,label='LT: {:.1f}'.format(LTsup),zorder=1)
            ax[6].plot(rhosdown.index, rhosdown.rho400/1e-12,'k',
                       linewidth=2,label='LT: {:.1f}'.format(LTsdown),zorder=0)
            ax[6].legend(frameon=False,loc=0)
            ax[6].set_ylabel(r'$\mathsf{\rho,\ 400 km}$',fontsize=14)
            for k1 in range(7):
                ax[k1].set_xlim(adates.min(),adates.max()+pd.Timedelta('1D'))
                ax[k1].xaxis.set_ticks(adates)
                ax[k1].grid(axis='x')
                ax[k1].yaxis.set_major_locator(MaxNLocator(nbins=4,prune='upper'))
            ax[-1].xaxis.set_ticklabels(adates.strftime('%y/%m/%d'),rotation=45)
            plt.subplots_adjust(hspace=0,top=0.95,right=0.95)
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
        density = density[abs(density.lat3)<=87]
        x = (density.index-pd.Timestamp('2004-12-31'))/pd.Timedelta('1D')
        y = density.lat
        z = density.rho400

        fig, ax = plt.subplots(3,2,sharex=True, sharey=True)
        plt.sca(ax[0,0])
        plt.tricontourf(x,y,z,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.title('plt.tricontourf')
        plt.ylabel('Lat')

        # dense grids
        xg, yg = np.meshgrid(np.arange(0,10,0.021),np.arange(-90,90,0.1))
        plt.sca(ax[1,0])
        zg = griddata((x*10,y),z,(xg*10,yg),method='linear',rescale=False)
        plt.contourf(xg,yg,zg,levels=np.arange(0,8,0.5)*1e-12)
        plt.xlim(1,10)
        plt.ylim(-90,90)
        plt.yticks(np.arange(-90,91,45))
        plt.ylabel('$\Delta$t=0.021D, $\Delta$lat=0.1')
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
        if False: # test champlt and gracelt
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
        if True: #  data preparation
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
                        smindex = get_index(k2+pd.TimedeltaIndex(np.arange(-lrday,lrday),'D'))
                        smindex = smindex.reindex(rhot.index,method='ffill')
                        rhot = rhot[smindex.Kp<=40]
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
        """ superposed epoch analysis
        """
        import gc
        lrday = 5
        sblist = get_sblist()
        sblist = sblist['2001-1-1':'2011-1-1']
        data = [pd.DataFrame(),pd.DataFrame()]
        for k1, k2 in enumerate(['away-toward', 'toward-away']):
            sbtmp = sblist[sblist.sbtype==k2]
            for idate in sbtmp.index:
                datelist = idate + pd.TimedeltaIndex(np.arange(10)-lrday, 'd')
                for satellite in ['champ','grace']:
                    density = get_density_dates(datelist,satellite=satellite)
                    if density.empty:
                        continue
                    density = density.add_updown()
                    for updown in ['up','down']:
                        d = density.orbit_mean(updown=updown)
                        d['epochday'] = d.index - idate
                        data[k1] = data[k1].append(d)
        pd.to_pickle(data,'/data/tmp/t6.dat')
        gc.collect()
    def f14():
        return


#--------------------------#
    #f5()
    #f6()
    #f7()
    #f8()
    #f9()
    #f10()
    #get_lt()
    f13()
    import gc
    gc.collect()
