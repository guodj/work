#-------------------------------------------------------------------------------
#
# By Dongjie, USTC, on Mon Sep 19 21:54:35 CST 2016
#
# Class for the GOCE densities and winds
#
# containe:
#
#       class: GoceData,
#           LT_median: Calculate the median local time of the ascending and
#                   descending orbits.
#           orbit_mean: Calculate the orbit mean longitude, altitude, local
#                   time, density.
#           satellite_position_lt_lat: show the satellite location in LT-LAT
#                   or MLT-MLAT coordinates
#           contourf_date_lat: Contourf of rho or winds as a function of
#                   date and latitude
#           polar_quiver_wind: wind vectors in the polar lat-lon coordinates.
#-------------------------------------------------------------------------------
# Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import myfunctions as mf
from scipy.interpolate import griddata

GOCEDATADIR = os.environ.get('DATAPATH') + 'GOCE/'
class GoceData(pd.DataFrame):
    """Class to open, manipulate and visualize the GOCE data files
    """

    def __init__(self, btime, etime, *args, **kwargs):
        super(GoceData, self).__init__(*args, **kwargs)
        self._read(btime, etime)
        if not self.empty:
            self.variables = tuple(self.columns)
            self.daterange = (self.index[0].strftime('%Y-%m-%d %H:%M:%S'),
                              self.index[-1].strftime('%Y-%m-%d %H:%M:%S'))

    def _read(self, bdate,edate):
        """ get goce data during specified dates.

        Args:
            bdate, edate: string or pd.Timestamp
        Returns:
            goce density indexed with datetime. the columns are:
            alt, long, lat, LT, arglat, rho, cr_wnd_e, cr_wnd_n, cr_wnd_u
        """
        global GOCEDATADIR
        bdate = pd.Timestamp(bdate)
        edate = pd.Timestamp(edate)
        dates = pd.date_range(bdate.date(),edate.date()+pd.Timedelta('1D'))
        yearmonth = np.unique(dates.strftime('%Y-%m'))
        yearmonth =pd.to_datetime(yearmonth)
        fname = [GOCEDATADIR+'/goce_denswind_v1_3_{:s}-{:s}.txt'.format(
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
            goce_data = goce_data[bdate:edate]
            for k0 in goce_data:
                self[k0] = goce_data[k0]
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
        fp = np.array((self.lat>=-30) &(self.lat<=30))
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
        # ascending or descending orbit?
        tmp = self[isup].copy() if updown =='up' else self[isdown].copy()
        tmp = tmp[(tmp.lat>=lats[0]) & (tmp.lat<=lats[1])]  #  which latitudes?
        tmp['float_time'] = (
                tmp.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
        # Below is a good method to calculate orbit number
        # Think carefully about the time interval between two orbits
        tmp['difft'] = np.insert(
                np.diff(tmp.index)/pd.Timedelta('1h'), 0, np.nan)
        tmp['orbitn'] = 0
        tmp.loc[tmp.difft>0.45,'orbitn']=1
        tmp['orbitn'] = tmp['orbitn'].cumsum()

        # There are at least half of the expected points
        maxn = tmp.orbitn.max()
        dp = len(tmp)
        grouped = tmp.groupby('orbitn').filter(
                lambda x: len(x.index)>dp/maxn/2).groupby('orbitn')
        result = grouped.agg(
                {'float_time': np.nanmean,
                 'long':np.nanmedian,
                 'alt':np.nanmean,
                 'LT':np.nanmedian,
                 'rho':np.nanmean})
        result['float_time'] = (result['float_time']*24*3600).round()
        result['datetime'] = (pd.Timestamp('2000-1-1') +
                              pd.TimedeltaIndex(result.float_time,'S'))
        result = result.set_index('datetime')
        result = result.drop('float_time',axis=1)
        return result


    def satellite_position_lt_lat(self, mag=False, ns='N'):
        """ Show the lt and lat positions of the satellite in a polar
        coordinate.

        Input:
            mag: if True, for MLT and Mlat position
            ns: N or S for North and South hemispheres, respectively

        Output:
            hcup, hcdown: scatter handles for the up and down orbits,
                respectively.
        """
        if self.empty:
            return
        lt = self['LT']
        lat = self['lat']
        if mag:
            from apexpy import Apex
            gm = Apex()
            mlat,mlt = gm.convert(
                    self['lat'], self['long'], 'geo', 'mlt', datetime=self.index)
            lt, lat = mlt, mlat
        ct = lat>0 if ns is 'N' else lat<0
        theta = lt[ct]/12*np.pi
        r = 90 - abs(lat[ct])
        hc = plt.scatter(theta, r, linewidths=0)
        return hc


    def contourf_date_lat(self, ax, whichcolumn='rho', **kwargs):
        """ A contourf of multiple-day density versus date and latitude.

        Args:
            ax: axis handle
            whichcolumn: string, 'rho', 'cr_wnd_e', 'cr_wnd_n','cr_wnd_u'.
            **kwargs: for contourf
        Return:
            hc: handle of the contourf plot
        """
        from matplotlib.ticker import AutoMinorLocator
        from scipy.interpolate import griddata
        import matplotlib.dates as mdates
        from matplotlib.ticker import AutoMinorLocator
        if not self.empty:
            self['epochday'] = (self.index-pd.Timestamp('2000-1-1'))/pd.Timedelta('1D')
            btime = self['epochday'].min()
            etime = self['epochday'].max()

            ut0 = np.arange(np.floor(btime), np.floor(etime)+1+0.1/24, 0.1/24)
            lat0 = np.arange(0,361,1)
            ut, lat = np.meshgrid(ut0, lat0)
            rho = griddata((self['epochday'], self.arglat),
                           self[whichcolumn], (ut, lat),
                           method='linear', rescale=True)
            for index, k in enumerate(ut0):
                fp = abs(self['epochday']-k)<0.5/24
                if not fp.any():
                    rho[:,index]=np.nan

            hc = ax.contourf(ut, lat, rho, 10, **kwargs)

            ax.set_xlim(np.floor(btime),np.floor(etime)+1)
            ax.set_xticks(np.arange(np.floor(btime),np.floor(etime)+2))
            ax.set_xticklabels(pd.date_range(
                    self.index[0],
                    self.index[-1]+pd.Timedelta('1d')).
                    strftime('%j'))
            ax.set_ylim(0,360)
            ax.set_yticks(np.arange(0,361,90))
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_minor_locator(AutoMinorLocator(3))
            ax.tick_params(which='both', width=1.2)
            ax.tick_params(which='major', length=7)
            ax.tick_params(which='minor', length=4)
            ax.tick_params(which='both',direction='out')
            ax.set_xlabel('Day of Year: {:d}'
                          .format(self.index[0].year),fontsize=14)
            ax.set_ylabel('Argument of Latitude', fontsize=14)
            return hc

    def polar_quiver_wind(self, ax, ns='N'):
        '''
        Wind vector in lat-lon coordinates.
        For different map projections, xy winds should be different?
        '''
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
        winde = self.cr_wnd_e
        windn = self.cr_wnd_n
        # only appropriate for the npstere and spstere
        xwind = fc*winde*np.cos(lon/180*np.pi)-windn*np.sin(lon/180*np.pi)
        ywind = winde*np.sin(lon/180*np.pi)+fc*windn*np.cos(lon/180*np.pi)
        hq = m.quiver(np.array(lon),np.array(lat),xwind,ywind,color='blue',
                      scale=800, scale_units='inches',zorder=3, latlon=True)
        #plt.quiverkey(hq,1.05,1.05,100,'100 m/s',coordinates='axes',labelpos='E')
        #m.scatter(np.array(lon),np.array(lat),
        #          s=50, c=self.index.to_julian_date(),linewidths=0, zorder=4,latlon=True)
        return m

#END
#--------------------------------------------------------------------------------
#TEST
if __name__ == '__main__':
    # Test get_goce_data, satellite_position_lt_lat
    den = GoceData('2010-4-1','2010-4-3')
    print(den.variables)
    print(den.daterange)
    print(den.orbit_mean())
    den.LT_median()
    plt.figure()
    ax = plt.subplot(polar=True)
    hc = den.satellite_position_lt_lat(mag=False)
    # Set polar(lat, LT) coordinates
    ax.set_rmax(30)
    ax.set_rgrids(np.arange(10,31,10),['$80^\circ$','$70^\circ$','$60^\circ$'],fontsize=14)
    ax.set_theta_zero_location('S')
    ax.set_thetagrids(np.arange(0,361,90),[0,6,12,18],fontsize=14,frac=1.05)
    # Test contourf_date_lat
    plt.figure()
    ax = plt.subplot()
    den.contourf_date_lat(ax,whichcolumn='cr_wnd_n')
    plt.tight_layout()
    plt.show()
    # Test polar_quiver_wind()
    plt.figure()
    ax = plt.subplot()
    wind = GoceData('2010-4-1','2010-4-1 1:30:0')
    wind.polar_quiver_wind(ax, ns='S')
    plt.show()
    #--------------------------------------------------
    # case study
    #    wind = get_goce_data('2010-5-29 6:00:00','2010-5-29 18:0:0')
    #    diffarglat = np.insert(np.diff(wind.arglat), 0, 0)
    #    wind['orbitn'] = 0
    #    wind.loc[diffarglat<0, 'orbitn']=1
    #    wind['orbitn'] = wind['orbitn'].cumsum()
    #    nm = wind.orbitn.max()+1
    #    fig,ax = plt.subplots(2,nm)
    #    for k0 in range(nm):
    #        plt.sca(ax[0, k0])
    #        tmp = GoceData(wind[wind.orbitn==k0])
    #        tmp.polar_quiver_wind(ax,ns='N')
    #        plt.sca(ax[1, k0])
    #        tmp.polar_quiver_wind(ax,ns='S')
    #    plt.tight_layout()
    #    plt.show()
