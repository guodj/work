#!/home/gdj/anaconda3/bin/python3
#-*- coding: utf-8 -*-

__author__ = 'Dongjie Guo'


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def get_champ_wind(dates):
    """ Get CHAMP wind during *contiuous days*.

    Input:
    dates: contiuous days such as pd.date_range('2005-1-1','2005-1-31',freq='D')

    Output
    champ_wind: dataframe of champ wind
    """
    champ_wind = pd.DataFrame()

    dates = pd.to_datetime(dates)
    dates = dates[(dates>'2001-1-1') & (dates<'2011-1-1')]
    datesfile = pd.to_datetime(pd.unique(dates.strftime('%Y-%m')))

    fname = ['/data/wind_CHAMPGFZ/'
             'acceldrag_datarequest_{:s}-{:s}-01_31d.txt'.format(
                 k.strftime('%Y'),k.strftime('%m')) for k in datesfile]
    champ_wind_tmp = [pd.read_csv(
        fn,
        parse_dates=[[0,1]],
        index_col=[0],
        comment='#',
        delim_whitespace=True,
        header=None,
        names=['date','time','time_system',
               'altitude','longitude','latitude','LT',
               'east_wind','north_wind']) for fn in fname]
    if champ_wind_tmp:
        champ_wind = pd.concat(champ_wind_tmp)
        champ_wind = champ_wind[
                (champ_wind.index>=dates.min()) &
                (champ_wind.index<dates.max()+pd.Timedelta('1D'))]
    return champ_wind


def set_lt_lat(h, pole='N'):
    """change some default sets of a polar axis

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


def quiver_champ_wind(bdate, edate, pole='N', updown='both', color='b'):
    """ Draw a quiver plot on the current axis
        showing the CHAMP wind around the poles.
        It is not considered if 2 orbits exist.

    Input:
    bdate: beginning date
    edate: end date
    pole: north or south pole
    updown: select up, down or both kinds of orbits
    color: color of the quiver

    Output:
    hs: handle of the scatter representing the orbit.
    hq: handle of the quiver showing the champ wind.
    """
    bdate = pd.to_datetime(bdate)
    edate = pd.to_datetime(edate)
    dates = pd.date_range(bdate.date(), edate.date())
    champ_wind = get_champ_wind(dates)
    if not champ_wind.empty:
        champ_wind = champ_wind[bdate:edate]
        lat = champ_wind.latitude
        dlat = lat.diff()
        dlat.iloc[0] = dlat.iloc[1]
        if updown == 'up':
            champ_wind = champ_wind[dlat>0]
        elif updown == 'down':
            champ_wind = champ_wind[dlat<0]
        if pole == 'N':
            champ_wind = champ_wind[champ_wind.latitude>=0]
            theta = champ_wind.LT*np.pi/12
            radius = 90-champ_wind.latitude
            champ_wind['xwind'] = (-champ_wind.east_wind*np.sin(theta-np.pi/2)
                                   -champ_wind.north_wind*np.cos(theta-np.pi/2))
            champ_wind['ywind'] = (champ_wind.east_wind*np.cos(theta-np.pi/2)
                                   -champ_wind.north_wind*np.sin(theta-np.pi/2))
        elif pole =='S':
            champ_wind = champ_wind[champ_wind.latitude<=0]
            theta = champ_wind.LT*np.pi/12
            radius = 90+champ_wind.latitude
            champ_wind['xwind'] = (-champ_wind.east_wind*np.sin(theta-np.pi/2)
                                   +champ_wind.north_wind*np.cos(theta-np.pi/2))
            champ_wind['ywind'] = (champ_wind.east_wind*np.cos(theta-np.pi/2)
                                   +champ_wind.north_wind*np.sin(theta-np.pi/2))
        # -np.pi/2 because zero LT points southward
        hq = plt.quiver(theta,radius,champ_wind.xwind,champ_wind.ywind,
                 color=color,
                 scale_units='x',scale=200,# set arrow length
                 headwidth=0,headlength=0,headaxislength=0)# no arrow head
        hs = plt.scatter(theta,radius,s=5,c='k',linewidths=0)
        set_lt_lat(plt.gca(),pole=pole)
        return hs, hq
    return


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
    if rho: # for list, if it is empty, then it is false
        rho = pd.concat(rho)
        rho = rho.groupby(rho.index).first()
        return rho
    else:
        return pd.DataFrame()


def time_champ_orbit_updown(bdate, edate, updown='up'):
    """Get the opening and terminate times of an up or down orbits.

    Input:
    bdate: beginning time of the time series
    edate: terminate time of the time series
    updown: up or down orbits?

    Output:
    betime: opening and terminate times
    """
    bdate = pd.to_datetime(bdate)
    edate = pd.to_datetime(edate)
    dates = pd.date_range(bdate.date(), edate.date())
    champ_wind = get_champ_wind(dates)
    if not champ_wind.empty:
        champ_wind = champ_wind[bdate:edate]
        lat = champ_wind.latitude
        dlat = lat.diff()
        dlat.iloc[0] = dlat.iloc[1]
        # The above guarantees maxdifftime[0] is not at the beginning of
        # difftime. In fact, it even does not matter whether it is at the
        # beginning.
        if updown == 'up':
            champ_wind = champ_wind[dlat>0]
        elif updown == 'down':
            champ_wind = champ_wind[dlat<0]
        else:
            print('Wrong input of updown!!!')
            return
        difftime = pd.to_timedelta(np.diff(champ_wind.index))
        difftime = difftime.append(pd.to_timedelta(['10s']))
        maxdifftime = difftime>'0.5h'
        etime = champ_wind.index[maxdifftime]

        difftime = pd.to_timedelta(np.roll(difftime,1))
        maxdifftime = difftime>'0.5h'
        btime = champ_wind.index[maxdifftime]

        btime = btime.insert(0,champ_wind.index[0])
        # etime.insert is wrong
        etime = etime.append(pd.DatetimeIndex([champ_wind.index[-1]]))

        betime = pd.DataFrame([btime,etime],index=['btime','etime']).T
        betime = betime[(betime.etime-betime.btime<'1h') &
                        (betime.etime-betime.btime>'0.1h')]
        return betime
    return pd.DataFrame()


def get_index(dates):
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


def get_imf_AE_PC_5m(dates):
    dates = pd.to_datetime(dates)
    years = np.unique(dates.strftime('%Y'))
    fname = ['/data/IMF_AE_PC_01_12_5m/csv/'
             '{:s}_csv.lst'.format(y) for y in years]
    index = [pd.read_csv(
        fn,
        parse_dates=[0],
        index_col=[0]) for fn in fname if os.path.isfile(fn)]
    if index:
        index = pd.concat(index)
        fp = np.floor(index.index.to_julian_date()+0.5).isin(
            dates.to_julian_date()+0.5)
        index = index.loc[fp]
        index.loc[index.Bym==9999.99, 'Bym'] = np.nan
        index.loc[index.Bzm==9999.99, 'Bzm'] = np.nan
        index.loc[index.AE==99999, 'AE'] = np.nan
        index.loc[index.PC==999.99, 'PC'] = np.nan
        return index
    else:
        return pd.DataFrame()


if __name__ == '__main__':
    """champ wind
    """
    bdate ='2003-10-4 10:00:00'
    edate ='2003-10-4 16:00:00'
    betime = time_champ_orbit_updown(bdate, edate, updown='up')
    fig = plt.figure(figsize=[8,8.25])
    for k, color, orbit in zip([1,2,3,4],
                               ['k',[0.5,0,0.5],'r','b'],
                               ['orbit1','orbit2','orbit3','orbit4']):
        btime = betime.loc[k-1,'btime']
        etime = betime.loc[k-1,'etime']
        h = plt.subplot(2,2,k,polar=True)
        hs,hq = quiver_champ_wind(btime,etime,pole='N',updown='up',color=color)
        h.set_title('NH: {:s}'.format(orbit), va='bottom')
        h.set_rmax(60)
    fig = plt.figure(figsize=[8,8.25])
    for k, color, orbit in zip([1,2,3,4],
                               ['k',[0.5,0,0.5],'r','b'],
                               ['orbit1','orbit2','orbit3','orbit4']):
        btime = betime.loc[k-1,'btime']
        etime = betime.loc[k-1,'etime']
        h = plt.subplot(2,2,k,polar=True)
        hs,hq = quiver_champ_wind(btime,etime,pole='S',updown='up',color=color)
        h.set_title('SH: {:s}'.format(orbit), va='bottom')
        h.set_rmax(60)
    plt.show()


    """champ density
    """
    #bdate = pd.to_datetime('2003-10-4 10:00')
    #edate = pd.to_datetime('2003-10-4 16:00')
    #dates = pd.date_range(bdate.date(), edate.date())
    #champ_density = get_density_dates(dates, satellite='champ')
    #betime = time_champ_orbit_updown(bdate, edate, updown='up')
    #fig =plt.figure(figsize=[8.8,3.7])
    #for kk, ns in enumerate(['N', 'S']):
    #    ax = plt.subplot(1,2,kk+1)
    #    for k, color, orbit in zip([1,2,3,4],
    #                               ['k',[0.5,0,0.5],'r','b'],
    #                               ['orbit1','orbit2','orbit3','orbit4']):
    #        btime = betime.loc[k-1,'btime']
    #        etime = betime.loc[k-1,'etime']
    #        tmp_density =champ_density[btime:etime]
    #        hp = plt.plot(tmp_density.lat, tmp_density.rho400/1e-12, '-o',
    #                      color=color,linewidth=2,
    #                      markeredgewidth=0,label=orbit)
    #        if ns == 'N':
    #            plt.axis([30,90,1.6,3])
    #            plt.title('NH')
    #            plt.xlabel('Latitude',fontsize=14)
    #            plt.ylabel('CHAMP Density ($10^{-12} kg/m^3$)')

    #        elif ns == 'S':
    #            plt.axis([-90,-30,1.6,3])
    #            plt.title('SH')
    #            plt.xlabel('Latitude',fontsize=14)
    #    if ns == 'N':
    #        ax.legend(loc='best',frameon=False)
    #fig.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
    #plt.show()


    #    """ AE variation
    #    """
    #    ax = plt.subplot()
    #    bdate = pd.to_datetime('2003-10-4 10:00')
    #    edate = pd.to_datetime('2003-10-4 16:00')
    #    dates = pd.date_range(bdate.date(), edate.date())
    #    index = get_imf_AE_PC_5m(dates)
    #    index = index[bdate:edate]

    #    betime = time_champ_orbit_updown(bdate, edate, updown='up')
    #    for k, color, orbit in zip([1,2,3,4],
    #                               ['k',[0.5,0,0.5],'r','b'],
    #                               ['orbit1','orbit2','orbit3','orbit4']):
    #        btime = betime.loc[k-1,'btime']
    #        etime = betime.loc[k-1,'etime']
    #        tmp_index = index[btime:etime]
    #        tmp_time = (tmp_index.index.hour +
    #                    tmp_index.index.minute/60 +
    #                    tmp_index.index.second/3600)
    #        tmp_AE = tmp_index.AE
    #        plt.plot(tmp_time,tmp_AE,color=color,label=orbit,linewidth=4, alpha=1)
    #    ax.legend(frameon=False)

    #    tmp_time = index.index.hour+index.index.minute/60+index.index.second/3600
    #    tmp_AE = index.AE
    #    plt.plot(tmp_time, tmp_AE, 'k--',linewidth=2,zorder=0,dashes=[10,5,10,5])

    #    plt.axis([10,16,0,500])
    #    plt.xlabel('Hours of 2003-10-4',fontsize=14)
    #    plt.ylabel('AE',fontsize=14)
    #    plt.grid()
    #    plt.show()
