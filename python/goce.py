#--------------------------------------------------------------------------------
# v 1.0
#
# By Dongjie, USTC, on Mon Sep 19 21:54:35 CST 2016
#
# Class for the GOCE densities.
# Also include a function to read data from file
#
# containe:
#       get_goce_data: Get density from GOCE satellites
#
#       class: GoceDensity,
#           print_variable_name: Print the column names
#
#           print_dates: Print the data dates
#
#           LT_median: Calculate the median local time of the ascending and
#               descending orbits.
#
#           add_updown: Add 2 columns that show the ascending and decending orbits
#
#           orbit_mean: Calculate the orbit mean longitude, altitude, local time,
#                   density.
#
#           satellite_position_lt_lat: show the satellite location in LT-LAT
#                   or MLT-MLAT coordinates
#
#           contourf_date_lat: Contourf of rho or winds as a function of
#                   date and latitude
#
#       ..........
#--------------------------------------------------------------------------------
# Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata

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
        return GoceDensity(goce_data)
    else:
        return GoceDensity()


class GoceDensity(pd.DataFrame):

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
        rho = rho[(rho.lat>=-30) &(rho.lat<=30)]
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


    def add_updown(self):
        """ Add 'isup' and 'isdown' columns to self
        Note that the function is appropriate for continuous data

        Returns:
            self added with columns 'isup' and 'isdown'

        Note:
            The results may be inappropriate near poles.
            Some bugs may exist at the data gap.
        """
        if not self.empty:
            lat = self.lat
            dlat = lat.diff()
            dlat.iloc[0] = dlat.iloc[1]
            self['isup'] = (dlat > 0)
            self['isdown'] = (dlat < 0)
            self = self[(self.isup) | (self.isdown)]
            return GoceDensity(self)


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
        tmp = self[self.isup] if updown =='up' else self[self.isdown]  # ascending or descending orbit?
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
            mag: if True, for MLT and Mlat position, Note that
                False is not supported now
            ns: N or S for North and South hemispheres, respectively

        Output:
            hcup, hcdown: scatter handles for the up and down orbits,
                respectively.
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


    def contourf_date_lat(self, ax, whichcolumn='rho',
                          updown='up', **kwargs):
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
            self['epochday'] = (self.index-self.index.min())/pd.Timedelta('1D')
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
                    strftime('%m-%d'),rotation=45)
            ax.set_ylim(0,360)
            ax.set_yticks(np.arange(0,361,90))
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_minor_locator(AutoMinorLocator(3))
            ax.tick_params(which='both', width=1.2)
            ax.tick_params(which='major', length=7)
            ax.tick_params(which='minor', length=4)
            ax.tick_params(which='both',direction='out')
            ax.set_xlabel('Date of Year: {:d}'
                          .format(self.index[0].year),fontsize=14)
            ax.set_ylabel('Argument of Latitude', fontsize=14)
            return hc


#END
#--------------------------------------------------------------------------------
#TEST
if __name__ == '__main__':
    den = get_goce_data(pd.date_range('2010-4-1','2010-4-30'))
    den.print_variable_name()
    den.print_dates()
    ax = plt.subplot(polar=True)
    hc = den.satellite_position_lt_lat(mag=False)
    #--------------------------------------------------------------------------------
    # Set polar(lat, LT) coordinates
    ax.set_rmax(30)
    ax.set_rgrids(np.arange(10,31,10),['$80^\circ$','$70^\circ$','$60^\circ$'],fontsize=14)
    ax.set_theta_zero_location('S')
    ax.set_thetagrids(np.arange(0,361,90),[0,6,12,18],fontsize=14,frac=1.05)
    #--------------------------------------------------------------------------------
    ax = plt.subplot()
    den.contourf_date_lat(ax,whichcolumn='cr_wnd_e')
    plt.tight_layout()
    plt.show()

