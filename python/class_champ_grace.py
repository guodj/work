"""
Class for the champ and grace density
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import griddata

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
        return output


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


#------------------------end class-------------------------------------#
