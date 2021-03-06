#-------------------------------------------------------------------------------
#
# By Dongjie, USTC/UM, on Fri Sep 16 23:23:39 CST 2016
#
# Class for the CHAMP and GRACE 3-degree densities and winds.
# Also include functions to read data from file
#
# Contain:
#
#       class: ChampDensity,
#           LT_median: Calculate the median local time of the ascending and
#                   descending orbits.
#
#           orbit_mean: Calculate the orbit mean longitude, height, local time,
#                   rho, rho400.
#
#           satellite_position_lt_lat: show the satellite location in LT-LAT
#                   coordinates
#
#           contourf_date_lat: Contourf of rho, rho400... as a function of date
#                   and lat
# Change:
#        1, Include ChampWind class
#        2, Use __init__, remove get_champ_grace_density, get_champ_wind
#        3, Add attributes variables and daterange, remove methods
#           print_variable_name and print_dates
#
#-------------------------------------------------------------------------------

# Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import myfunctions as mf
import os

DATADIR = os.environ.get('DATAPATH')
class ChampDensity(pd.DataFrame):
    """
    Class to open, manipulate and visualize the CHAMP or/and GRACE
    density files.
    """

    def __init__(self, btime, etime, satellite='champ', variables=-1,
                 *args, **kwargs):
        super(ChampDensity, self).__init__(*args, **kwargs)
        self._read(btime, etime, satellite, variables)
        # if not self.empty:
        #     self.variables = tuple(self.columns)

    def _read(self, bdate, edate, satellite, variables):
        """ get champ or grace density data during specified period.

        Args:
            bdate, edate: string or pd.Timestamp
            satellite: 'champ' or 'grace'
            variables: subset of default.
        Returns:
            ChampDensity of champ or grace density indexed with datetime.
        """
        global DATADIR
        bdate = pd.Timestamp(bdate)
        edate = pd.Timestamp(edate)
        dates = pd.date_range(bdate.date(), edate.date(), freq='1D')
        # champ and grace data are in this date range
        dates = dates[(dates>'2001-1-1') & (dates<'2011-1-1')]
        column_names=[
                'year','doy','second',
                'lat3','lat','long','height','LT','Mlat','Mlong','MLT',
                'rho','rho400','rho410','msis_rho',
                'uncertainty','data_points','points_required',
                'drag_coefficient']
        if variables == -1:
            variables = ['lat3', 'lat', 'long', 'height', 'LT',
                         'Mlat', 'Mlong', 'MLT',
                         'rho', 'rho400', 'rho410', 'msis_rho']
        variables = list(variables)
        variables.append('second')
        rho = []  # List to store data
        for k0 in dates:
            if satellite == 'champ':
                fname = (DATADIR + 'CHAMP/density/{:s}/ascii/'
                        'Density_3deg_{:s}_{:s}.ascii'.format(
                                k0.strftime('%Y'),
                                k0.strftime('%y'),
                                k0.strftime('%j')))
            if satellite == 'grace':
                fname = (DATADIR + 'Grace/{:s}/ascii/'
                        'Density_graceA_3deg_{:s}_{:s}.ascii'.format(
                                k0.strftime('%Y'),
                                k0.strftime('%y'),
                                k0.strftime('%j')))
            if os.path.isfile(fname):
                tmp = pd.read_csv(
                        fname,
                        delim_whitespace=True,
                        skiprows=2,
                        header=None,
                        usecols=variables,
                        names=column_names)
                if not tmp.empty:
                    tmp['date'] = k0 + pd.TimedeltaIndex(tmp.second, 's')
                    tmp.set_index('date', drop=True, append=False, inplace=True)
                    tmp.drop('second', axis=1, inplace=True)
                    rho.append(tmp)
        if rho:
            rho = pd.concat(rho)
            # Exclude duplicate points
            # pd.DataFrame.drop_duplicates() has something wrong
            rho = rho.groupby(rho.index).first()
            #rho.drop_duplicates(inplace=True)
            rho = rho[bdate:edate]
            for k0 in rho:
                self[k0] = rho[k0]
        return

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
        isup, isdown = mf.updown(self.lat)
        fp = np.array((self.lat3>=-30) &(self.lat3<=30))
        rho, isup, isdown = self[fp].copy(), isup[fp], isdown[fp]
        if rho.empty:
            return output
        grouped = rho.groupby(isup)['LT']
        for name, group in grouped:
            k0 = 0 if name is True else 1
            group1 = group
            if group1.max()-group1.min()>22:
                group1[group1<2] = group1[group1<2]+24
            output[k0] = np.median(group1)
            output[k0] = output[k0]%24
        print('Ascending LT: %4.1f, Descending LT: %4.1f'%(output[0],output[1]))
        return output

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
        isup, isdown = mf.updown(self.lat)
        tmp = self[isup] if updown =='up' else self[isdown]  # ascending or descending orbit?
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


    def satellite_position_lt_lat(self, ax, nlat=90, slat=0,
                                  mag=False, plot_type='polar',
                                  *args, **kwargs):
        """ Show the (M)lt and (M)lat positions of the satellite in a polar
        or rectangular coordinate axis.

        Input:
            ax: which axis to draw orbit on
            nlat: northward latitude limit (default 90)
            slat: southward latitude limit (default 0)
            mag: if True, for MLT and Mlat position
        Return:
            hc: handle of the scatter
        """
        if self.empty:
            return
        lt='MLT' if mag else 'LT'
        lat='Mlat' if mag else 'lat'
        ct = (self[lat]>slat) & (self[lat]<nlat)
        if 'pol' in plot_type.lower():
            if nlat*slat<0:
                print('Error: For polar plot, nlat and slat'
                      ' should have the same signs')
            csign = 1 if nlat>0 else -1
            theta = self.loc[ct,lt]/12*np.pi
            r = 90-csign*self.loc[ct,lat]
            hc = ax.scatter(theta, r, linewidths=0, *args, **kwargs)
        if 'rec' in plot_type.lower():
            hc = ax.scatter(self.loc[ct, lt], self.loc[ct, lat], linewidths=0,
                            *args, **kwargs)
        return hc


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
        from scipy.interpolate import griddata
        if not self.empty:
            self['epochday'] = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
            btime = self['epochday'].min()
            etime = self['epochday'].max()

            isup, isdown = mf.updown(self.lat)
            tmp = self[isup] if updown is 'up' else self[isdown]

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


class ChampWind(ChampDensity):
    """Class to open, manipulate and visualize the CHAMP
    wind files.
    """
    def __init__(self, btime, etime, variables=-1, *args, **kwargs):
        super(ChampDensity, self).__init__(*args, **kwargs)
        self._read(btime, etime, variables)
        self.variables = tuple(self.columns)

    def _read(self, bdate, edate, variables):
        # Get champ winds during 'bdate' and 'edate'
        # variables is a list(tuple)
        global DATADIR
        bdate = pd.Timestamp(bdate)
        edate = pd.Timestamp(edate)
        dates = pd.date_range(bdate.date(),edate.date(),freq='1D')
        # champ and grace data are in this date range
        dates = dates[(dates>'2001-1-1') & (dates<'2011-1-1')]
        column_names=[
                'year','doy','second',
                'lat3','lat','long','height','LT',
                'wind', 'winde', 'windn','t1','t2','t3']
        if variables == -1:
            variables = ['lat3', 'lat', 'long', 'height', 'LT',
                         'wind', 'winde', 'windn']
        variables = list(variables)
        variables.append('second')
        wind = []  # List to store data
        for k0 in dates:
            fname = (DATADIR + 'CHAMP/winds/{:s}/'
                    'Wind_3deg_{:s}_{:s}.ascii'.format(
                            k0.strftime('%Y'),
                            k0.strftime('%y'),
                            k0.strftime('%j')))
            if os.path.isfile(fname):
                tmp = pd.read_csv(
                        fname,
                        delim_whitespace=True,
                        skiprows=2,
                        header=None,
                        usecols=variables,
                        names=column_names)
                if not tmp.empty:
                    tmp['date'] = k0 + pd.TimedeltaIndex(tmp.second, 's')
                    tmp.set_index('date', drop=True, append=False, inplace=True)
                    tmp.drop('second', axis=1, inplace=True)
                    wind.append(tmp)
        if wind:
            wind = pd.concat(wind)
            # Exclude duplicate points
            # pd.DataFrame.drop_duplicates() has something wrong
            wind = wind.groupby(wind.index).first()
            wind = wind[bdate:edate]
            for k0 in wind:
                self[k0] = wind[k0]
        return

    def satellite_position_lt_lat(self, ax, nlat=90, slat=0,
                                  mag=False, plot_type='polar',
                                  *args, **kwargs):
        """ Show the (M)lt and (M)lat positions of the satellite in a polar
        or rectangular coordinate axis.

        Input:
            ax: which axis to draw orbit on
            nlat: northward latitude limit (default 90)
            slat: southward latitude limit (default 0)
            mag: if True, for MLT and Mlat position
        Return:
            hc: handle of the scatter
        """
        if self.empty:
            return
        tmp = self
        if mag:
            from apexpy import Apex
            import datetime as dt
            a = Apex(date=self.index.year.mean())
            mlat, mlt = a.convert(tmp.lat, tmp.long, 'geo','mlt',
                                  height=tmp.height, datetime=tmp.index)
            tmp['MLT'] = mlt
            tmp['Mlat'] = mlat
        ltp='MLT' if mag else 'LT'
        latp='Mlat' if mag else 'lat'
        ct = (self[latp]>slat) & (self[latp]<nlat)
        if 'pol' in plot_type.lower():
            if nlat*slat<0:
                print('Error: For polar plot, nlat and slat'
                      ' should have the same signs')
            csign = 1 if nlat>0 else -1
            theta = tmp.loc[ct,ltp]/12*np.pi
            r = 90-csign*tmp.loc[ct,latp]
            hc = ax.scatter(theta, r, linewidths=0, *args, **kwargs)
        if 'rec' in plot_type.lower():
            hc = ax.scatter(tmp.loc[ct, ltp], tmp.loc[ct, latp], linewidths=0,
                            *args, **kwargs)
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
        from scipy.interpolate import griddata
        if not self.empty:
            #self['epochday'] = (self.index-self.index.min())/pd.Timedelta('1D')
            self['epochday'] = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
            btime = self['epochday'].min()
            etime = self['epochday'].max()

            isup, isdown = mf.updown(self.lat)
            tmp = self[isup] if updown is 'up' else self[isdown]

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
#-------------------------------------------------------------------------------
# for test
if __name__=='__main__':
    #------------------------------------------------------------
    # Check whether satellite_position_lt_lat results
    # from ChampWind and ChampDensity are the same.
    wind = ChampWind('2006-1-1 12:00:00','2006-1-1 13:30:00')
    ax = plt.subplot()
    wind.satellite_position_lt_lat(
            ax, mag=True, nlat=90, slat=-90, plot_type='rec', color='b')
    density = ChampDensity('2006-1-1 12:00:00','2006-1-1 13:30:00')
    density.satellite_position_lt_lat(
            ax, mag=True, nlat=90, slat=-90, plot_type='rec', color='r')
    plt.show()
    #------------------------------------------------------------
    # Test ChampWind.contourf_date_lat.
    #    wind = ChampWind('2005-1-1','2005-1-10 23:59:59')
    #    ax = plt.subplot()
    #    wind.contourf_date_lat(ax,whichcolumn='wind')
    #    plt.show()
    #------------------------------------------------------------
    # Test ChampWind.polar_quiver_wind.
    #    wind = ChampWind('2005-11-1','2005-11-1 5:0:0')
    #    wind['arglat'] = mf.lat2arglat(wind.lat)
    #    diffarglat = np.insert(np.diff(wind.arglat), 0, 0)
    #    wind['orbitn'] = 0
    #    wind.loc[diffarglat<0, 'orbitn']=1
    #    wind['orbitn'] = wind['orbitn'].cumsum()
    #    nm = wind.orbitn.max()+1
    #    fig,ax = plt.subplots(nm,2)
    #    for k0 in range(nm):
    #        plt.sca(ax[k0,0])
    #        tmp = wind[wind.orbitn==k0]
    #        tmp.__class__ = ChampWind
    #        tmp.polar_quiver_wind(ax,ns='N')
    #        plt.sca(ax[k0,1])
    #        tmp.polar_quiver_wind(ax,ns='S')
    #    plt.tight_layout()
    #    plt.show()
    #------------------------------------------------------------
    # Test add_arglat
    #    wind = ChampWind('2005-11-11','2005-11-14')
    #    wind['arglat'] = mf.lat2arglat(wind.lat)
    #    plt.plot(wind.index,wind.arglat)
    #    plt.show()
    # Test attributes variables and daterange
    #    a = ChampDensity('2005-1-1', '2005-1-2')
    #    print(a.variables)
    #    print(a.daterange)
    # Test satellite_position_lt_lat
    #    a = ChampDensity('2005-1-1', '2005-1-2')
    #    ax = plt.subplot()
    #    a.satellite_position_lt_lat(ax, nlat=90, slat=-90, plot_type='rec')
    #    plt.show()

