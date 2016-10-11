#--------------------------------------------------------------------------------
# v 1.0
#
# By Dongjie, USTC, on Fri Sep 16 23:23:39 CST 2016
#
# Class for the CHAMP and GRACE 3-degree densities and winds.
# Also include functions to read data from file
#
# Contain:
#       get_champ_grace_data: Get density from CHAMP or GRACE satellites
#
#       class: ChampDensity,
#           print_variable_name: Print the column names
#
#           print_dates: Print the data dates
#
#           LT_median: Calculate the median local time of the ascending and
#               descending orbits.
#
#           add_updown: Add 2 columns that show the ascending and decending orbits
#
#           orbit_mean: Calculate the orbit mean longitude, height, local time,
#                   rho, rho400.
#
#           satellite_position_lt_lat: show the satellite location in LT-LAT
#                   coordinates
#
#           data_gap_nan: Add nan to data gap, used for plot. There is a data gap
#                   when the time interval is greater than 2 hours
#
#           contourf_date_lat: Contourf of rho, rho400... as a function of date and
#                   lat
# Change:
#        Include wind handling, on Sat Sep 24 02:18:03 CST 2016
#            1, get_champ_wind, Get Champ wind
#            2, class ChampWind, Subclass of ChampDensity in order to use its function.
#
#       ..........
#--------------------------------------------------------------------------------

# Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata

DATADIR = '/home/guod/data/'
def get_champ_grace_density(
        bdate,edate,satellite='champ', variables=(
            'lat3', 'lat', 'long', 'height', 'LT',
            'Mlat', 'Mlong', 'MLT',
            'rho', 'rho400', 'rho410', 'msis_rho')):
    """ get champ or grace density data during specified dates.

    Args:
        bdate, edate: string or pd.Timestamp
        satellite: 'champ' or 'grace'
    Returns:
        dataframe of champ or grace density indexed with datetime. the columns
        are:lat3, lat, long, height, LT, Mlat, Mlong, MLT, rho, rho400,rho410,
        msis_rho, ...
    """
    global DATADIR
    bdate = pd.Timestamp(bdate)
    edate = pd.Timestamp(edate)
    dates = pd.date_range(bdate.date(),edate.date()+pd.Timedelta('1D'))
    dates = dates[(dates>'2001-1-1') & (dates<'2011-1-1')]
    variables = list(variables)
    variables.append('date')
    # champ and grace data are in this date range
    if satellite == 'champ':
        fname = [DATADIR+'CHAMP/density/csv/{:s}/'
                 'Density_3deg_{:s}_{:s}.ascii'.format(
                     k.strftime('%Y'),
                     k.strftime('%y'),
                     k.strftime('%j')) for k in dates]
    elif satellite == 'grace':
        fname = [DATADIR+'Grace/csv/{:s}/ascii/'
                 'Density_graceA_3deg_{:s}_{:s}.ascii'.format(
                     k.strftime('%Y'),
                     k.strftime('%y'),
                     k.strftime('%j')) for k in dates]
    rho = [pd.read_csv(
        fn,
        usecols = variables,
        parse_dates=['date'],
        index_col=['date']) for fn in fname if os.path.isfile(fn)]
    if rho:
        rho = pd.concat(rho)
        # Exclude duplicate points
        # pd.DataFrame.drop_duplicates() has something wrong
        rho = rho.groupby(rho.index).first()
        rho = rho[bdate:edate]
        return ChampDensity(rho)
    else:
        return ChampDensity()


def get_champ_wind(
        bdate, edate, variables=('lat3','lat','long','height','LT','wind','winde','windn')):
    # Get champ winds during 'bdate' and 'edate'
    # variables is a list(tuple)
    global DATADIR
    bdate = pd.Timestamp(bdate)
    edate = pd.Timestamp(edate)
    dates = pd.date_range(bdate.date(),edate.date(),freq='1D')
    fname = [DATADIR+'CHAMP/winds/csv/{:s}/Wind_3deg_{:s}_{:s}.ascii'.
            format(k.strftime('%Y'), k.strftime('%y'), k.strftime('%j'))
            for k in dates]
    variables = list(variables)
    variables.extend(['date'])
    wind = [pd.read_csv(fn,parse_dates=['date'],index_col=['date'],
                        usecols=variables,squeeze=True)
            for fn in fname if os.path.isfile(fn)]
    if wind:
        wind = pd.concat(wind)
        wind = wind[bdate:edate]
        return ChampWind(wind)
    else:
        return ChampWind()

class ChampDensity(pd.DataFrame):

    def print_variable_name(self):
        # Print column names in self
        if self.empty:
            return
        for k00,k0 in enumerate(self.columns):
            print(k00,': ',k0)


    def print_dates(self):
        # Print dates of the data
        if self.empty:
            return
        print(pd.DatetimeIndex(np.unique(self.index.date)))


    def LT_median(self):
        """ Get the local time of the ascending and descending satellite orbit.
        Only for a short-period continuous period.

        Returns:
            upLT, downLT: Local times for ascending and descending orbits.
        Note:
            Only low latitudes (-30<lat<30) are considered.
        """
        output = [np.nan,np.nan]
        if self.empty:
            return output
        rho = self.add_updown()
        rho = rho[(rho.lat3>=-30) &(rho.lat3<=30)]
        if rho.empty:
            return output
        grouped = rho.groupby(rho.isup)['LT']
        for name, group in grouped:
            k0 = 0 if name is True else 1
            group1 = group
            if group1.max()-group1.min()>22:
                group1[group1<2] = group1[group1<2]+24
            output[k0] = np.median(group1)
            output[k0] = output[k0]%24
        print('Ascending LT: %4.1f, Descending LT: %4.1f'%(output[0],output[1]))
        return output

    def add_updown(self, whichlat='lat3'):
        """ Add 'isup' and 'isdown' columns to self
        Note that the function is appropriate for continuous data

        Args:
            whichlat: data column name used to seperate up and down orbits
        Returns:
            self added with columns 'isup' and 'isdown'

        Note:
            The results may be inappropriate near poles.
            Some bugs may exist at the data gap.
            --------------------
            use a = a.add_updown(), not a.add_updown()
            --------------------
            Use updown() in myfunctions.py, that may be better
        """
        if not self.empty:
            lat = self[whichlat]
            dlat = lat.diff()
            dlat.iloc[0] = dlat.iloc[1]
            if whichlat=='lat3':
                fp = (np.abs(dlat)!=3) & (np.abs(dlat)!=0)
                fpid = np.argwhere(fp).reshape(-1)[:-1]
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

    def add_arglat(self):
        # This algorithm is wrong. But let it be.
        # Argument of latitude is the angle from the ascending node
        # to the satellite position along the orbital path.
        if not 'isup' in self:
            self = self.add_updown()
        self['arglat'] = self.lat3.copy()
        self.loc[self.isup,'arglat'] = (360+self.loc[self.isup,'arglat'])%360
        self.loc[self.isdown,'arglat'] = 180-self.loc[self.isdown,'arglat']
        return ChampDensity(self)

    def orbit_mean(self,lats=(-85,85),updown='up'):
        """ Get the orbit mean density during specified latitude range and
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

        # There are at least half of the expected points
        maxn = tmp.orbitn.max()
        dp = len(tmp)
        grouped = tmp.groupby('orbitn').filter(
                lambda x: len(x.index)>dp/maxn/2).groupby('orbitn')
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


    def satellite_position_lt_lat(self, mag=False, ns='N'):
        """ Show the (M)lt and (M)lat positions of the satellite in a polar
        coordinate.

        Input:
            mag: if True, for MLT and Mlat position
            ns: N or S for North and South hemispheres, respectively
        """
        if self.empty:
            return
        lt='MLT' if mag else 'LT'
        lat='Mlat' if mag else 'lat'
        ct = self[lat]>0 if ns is 'N' else self[lat]<0
        theta = self.loc[ct,lt]/12*np.pi
        r = 90 - abs(self.loc[ct,lat])
        hc = plt.scatter(theta, r, linewidths=0)
        return hc


    def data_gap_nan(self):
        """ If there are data gaps,insert nan, for plt.plot
        """
        tmp = self.index
        tmp1 = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1H')
        tmp2 = np.insert(np.diff(tmp1),0,0)
        tmp3 = np.argwhere(tmp2>2)
        tmp3 = tmp3.reshape((tmp3.size,))
        tmp = np.insert(tmp,tmp3,tmp[tmp3-1]+pd.Timedelta('0.5H'))
        return self.reindex(tmp)


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


    def get_lt_lat_density(self, whichdensity='rho400', latbin=10,
                           ltbin=8, interpolate=False, latn=180, ltn=48):
        """ Get lt, lat, data for polar contourf. This function is for a large
        data set that covers all the lt and lat ranges.

        Args:
            whichdensity: rho, rho400, rho410...
            latbin: lat window, harmonic of 90: 3,10,30....
            ltbin: local time window, harmonic of 24: 1,2,3,4,6,8...
            interpolate: bool, whether to interpolate data to higher resolution
            latn, ltn: for interpolation, bin numbers of lat and lt, integral
                multiple of 180 and 24 respectively
        Returns:
            lt, lat, density. All are 2D array

        Note: the sample number in each lat*lt grid is not considered
              maybe the interpolation is unnecessary
        """
        df = self[['lat', 'LT', whichdensity]]
        df['lat'] = pd.cut(df['lat'], bins=np.arange(-90,90+1e-10,latbin),
                     labels=np.arange(-90,90+1e-10,latbin)[:-1]+latbin/2,
                     include_lowest=True)
        df['LT'] = pd.cut(df['LT'], bins=np.arange(0,24+1e-10,ltbin),
                     labels=np.arange(0,24+1e-10,ltbin)[:-1]+ltbin/2,
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


    def get_lon_lat_density(self, whichdensity='rho400', latbin=10,
                            lonbin=45, interpolate=False,
                            latn=180, lonn=360):
        """ return lon, lat, density for polar contourf. This function is for a
        large data set that covers all the lon and lat ranges.

        Args:
            whichdensity: rho, rho400, rho410...
            latbin: lat window, harmonic of 90: 3,10...
            lonbin: lon window, harmonic of 360: 10,30...
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
                df['lat'], bins=np.arange(-90,90+1e-10,latbin),
                labels=np.arange(-90,90+1e-10,latbin)[:-1]+latbin/2,
                include_lowest=True)
        df['lon'] = pd.cut(
                df['lon'], bins=np.arange(-180,180+1e-10,lonbin),
                labels=np.arange(-180,180+1e-10,lonbin)[:-1]+lonbin/2,
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


class ChampWind(ChampDensity):
    def satellite_position_lt_lat(self, mag=False, ns='N'):
        """ Show the (M)lt and (M)lat positions of the satellite in a polar
        coordinate.

        Input:
            mag: if True, for MLT and Mlat position
            ns: N or S for North and South hemispheres, respectively
        """
        if self.empty:
            return
        tmp = self
        if mag:
            from apexpy import Apex
            import datetime as dt
            a = Apex(date=2005)
            mlat,mlt = a.convert(tmp.lat, tmp.long, 'geo','mlt', datetime=tmp.index)
            tmp['MLT'] = mlt
            tmp['Mlat'] = mlat
            #for k00,k0 in enumerate(tmp.index):
            #    mlat,mlt = a.convert(tmp.ix[k00,'lat'], tmp.ix[k00,'long'],
            #                         'geo','mlt', datetime=k0)
            #    tmp.ix[k00,'Mlat'] = mlat
            #    tmp.ix[k00,'MLT'] = mlt
        ltp='MLT' if mag else 'LT'
        latp='Mlat' if mag else 'lat'
        ct = tmp[latp]>0 if ns is 'N' else tmp[latp]<0
        theta = tmp.loc[ct,ltp]/12*np.pi
        r = 90 - abs(tmp.loc[ct,latp])
        hc = plt.scatter(theta, r, linewidths=0)
        return hc


    def contourf_date_lat(self, ax, whichcolumn='wind', updown='up', **kwargs):
        """ A contourf of multiple-day wind versus date and latitude.

        Args:
            ax: axis handle
            whichcolumn: string, 'wind', 'winde', 'windn'.
            updown: string, 'up' or 'down'
            **kwargs: for contourf
        Return:
            hc: handle of the contourf plot
        ----------------------------------------
        Note: x axis is days from '2000-1-1'
        """
        from matplotlib.ticker import AutoMinorLocator
        if not self.empty:
            #self['epochday'] = (self.index-self.index.min())/pd.Timedelta('1D')
            self['epochday'] = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
            btime = self['epochday'].min()
            etime = self['epochday'].max()

            self = self.add_updown()
            tmp = self[self.isup] if updown is 'up' else self[self.isdown]

            ut0 = np.arange(np.floor(btime), np.floor(etime)+1+0.5/24, 0.5/24)
            lat0 = np.arange(-90,91,3)
            ut, lat = np.meshgrid(ut0, lat0)
            windt = griddata((tmp['epochday'], tmp.lat), tmp[whichcolumn], (ut, lat),
                             method='linear', rescale=True)
            for index, k in enumerate(ut0):
                fp = abs(tmp['epochday']-k)<0.5/24
                if not fp.any():
                    windt[:,index]=np.nan
            hc = ax.contourf(
                    ut, lat, windt,
                    levels=np.linspace(np.nanpercentile(windt,1),np.nanpercentile(windt,99),11),
                    **kwargs)
            ax.set_xlim(np.floor(btime),np.floor(etime)+1)
            ax.set_xticks(np.arange(np.floor(btime),np.floor(etime)+2))
            ax.set_xticklabels(
                    pd.date_range(tmp.index[0],tmp.index[-1]+pd.Timedelta('1d')).strftime('%j'))
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
            return hc#, windt


    def polar_quiver_wind(self, ax, ns='N'):
        # Wind vector in lat-long coordinates.
        # For different map projections, the arithmetics to calculate xywind
        # are different
        if self.empty:
            return
        from mpl_toolkits.basemap import Basemap
        from apexpy import Apex
        # Creat polar coordinates
        projection,fc = ('npstere',1) if ns=='N' else ('spstere',-1)
        m = Basemap(projection=projection,boundinglat=fc*40,lon_0=0,resolution='l')
        m.drawcoastlines(color='gray',zorder=1)
        m.fillcontinents(color='lightgray',zorder=0)
        dt = self.index.min() + (self.index.max()-self.index.min())/2
        m.nightshade(dt,zorder=2)
        #m.drawparallels(np.arange(-80,81,20))
        #m.drawmeridians(np.arange(-180,181,60),labels=[1,1,1,1])

        # Calculate mlat and mlon
        lat_grid = np.arange(-90,91,10)
        lon_grid = np.arange(-180,181,10)
        lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)
        gm = Apex(date=2005)
        mlat,mlon = gm.convert(lat_grid,lon_grid,'geo','apex')
        hc1 = m.contour(lon_grid,lat_grid,mlat,levels=np.arange(-90,91,10),
                        colors='k', zorder=3, linestyles='dashed',
                        linewidths=1, latlon=True)
        # hc2 = m.contour(lon_grid,lat_grid,mlon,levels=np.arange(-180,181,45),
        #                 colors='k', zorder=3, linestyles='dashed', latlon=True)
        plt.clabel(hc1,inline=True,colors='k',fmt='%d')
        # plt.clabel(hc2,inline=True,colors='k',fmt='%d')

        # Calculate and plot x and y winds
        lat = self.lat
        lon = self.long
        wind = self.wind
        winde1 = self.winde
        winde = winde1*wind
        windn1 = self.windn
        windn = windn1*wind
        # only appropriate for the npstere and spstere
        xwind = fc*winde*np.cos(lon/180*np.pi)-windn*np.sin(lon/180*np.pi)
        ywind = winde*np.sin(lon/180*np.pi)+fc*windn*np.cos(lon/180*np.pi)
        hq = m.quiver(np.array(lon),np.array(lat),xwind,ywind,color='blue',
                      scale=300, scale_units='inches',zorder=3, latlon=True)
        #plt.quiverkey(hq,1.05,1.05,100,'100 m/s',coordinates='axes',labelpos='E')
        #m.scatter(np.array(lon),np.array(lat),
        #          s=50, c=self.index.to_julian_date(),linewidths=0, zorder=4,latlon=True)
        return m
# END
#--------------------------------------------------------------------------------
# for test
if __name__=='__main__':
    #------------------------------------------------------------
    #    # Test satellite_position_lt_lat.
    #    den = get_champ_grace_data(
    #           '2005-1-2','2005-1-3', variables=['lat','Mlat','LT','MLT','rho400','rho'])
    #    ax = plt.subplot(polar=True)
    #    hc = den.satellite_position_lt_lat(mag=True)
    #    #----------------------------------------
    #    # Set polar(lat, LT) coordinates
    #    ax.set_rmax(30)
    #    ax.set_rgrids(
    #           np.arange(10,31,10),['$80^\circ$','$70^\circ$','$60^\circ$'],fontsize=14)
    #    ax.set_theta_zero_location('S')
    #    ax.set_thetagrids(np.arange(0,361,90),[0,6,12,18],fontsize=14,frac=1.05)
    #    #----------------------------------------
    #    plt.show()
    #------------------------------------------------------------
    #    # Check whether satellite_position_lt_lat results
    #    # from ChampWind and ChampDensity are the same.
    #    wind = get_champ_wind('2005-1-1 12:00:00','2005-1-10 13:00:00')
    #    plt.figure()
    #    ax = plt.subplot(polar=True)
    #    wind.satellite_position_lt_lat(mag=True)
    #    density = get_champ_grace_data(pd.date_range('2005-1-1','2005-1-10'))
    #    plt.figure()
    #    ax = plt.subplot(polar=True)
    #    density.satellite_position_lt_lat(mag=True)
    #    plt.show()
    #------------------------------------------------------------
    #    # Test ChampWind.contourf_date_lat.
    #    wind = get_champ_wind('2005-1-1','2005-1-10')
    #    ax = plt.subplot()
    #    wind.contourf_date_lat(ax,whichcolumn='wind')
    #    plt.show()
    #------------------------------------------------------------
    # Test ChampWind.polar_quiver_wind.
    wind = get_champ_wind('2005-11-1','2005-11-1 3:0:0')
    wind = wind.add_arglat()
    diffarglat = np.insert(np.diff(wind.arglat), 0, 0)
    wind['orbitn'] = 0
    wind.loc[diffarglat<0, 'orbitn']=1
    wind['orbitn'] = wind['orbitn'].cumsum()
    nm = wind.orbitn.max()+1
    fig,ax = plt.subplots(nm,2)
    for k0 in range(nm):
        plt.sca(ax[k0,0])
        tmp = ChampWind(wind[wind.orbitn==k0])
        tmp.polar_quiver_wind(ax,ns='N')
        plt.sca(ax[k0,1])
        tmp.polar_quiver_wind(ax,ns='S')
    plt.tight_layout()
    plt.show()
    #------------------------------------------------------------
    # Test add_arglat
    #    wind = get_champ_wind('2005-11-11','2005-12-14')
    #    wind = wind.add_arglat()
    #    plt.plot(wind.index,wind.arglat)
    #    plt.show()
