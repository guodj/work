#--------------------------------------------------------------------------------
# Functions that may be used many times
# By Dongjie, start from Mon Oct 10 09:37:12 EDT 2016
#--------------------------------------------------------------------------------

# Global imports

import numpy as np
import pandas as pd

def updown(lat):
    '''
    Find ascending and descending orbits.
    'lat' is satellite continuous latitudes in the form of
    'pd.Series' with 'pd.DatetimeIndex' as index.
    Return np.ndarray with data type bool
    ----------------------------------------
    1, Data gaps are taken into consideration
    2, Inappropriate near poles, but this matters little.
    '''
    dlat = np.diff(lat)
    dlat = np.insert(dlat,0,dlat[0])
    # Set dlat at data gap as the next value.
    dt = (lat.index - pd.Timestamp('2000-1-1'))/pd.Timedelta('1m')
    ddt = np.diff(dt)
    ddt = np.insert(ddt,0,ddt[0])
    mddt = np.nanmedian(ddt)
    fp, = np.where(ddt > 10*mddt)
    if fp.size != 0:
        if fp.max() == ddt.size-1:
            fp = fp[:-1]
        dlat[fp] = dlat[fp+1]
    isup = dlat>0
    isdown = dlat<0
    return (isup, isdown)

def lat2arglat(lat):
    '''
    convert latitude to argument of latitude

    'lat' is satellite continuous latitudes in the form of
    'pd.Series' with 'pd.DatetimeIndex' as index.

    Note: only right for polar near circular orbit
    '''
    arglat = np.array(lat.copy()*np.nan)
    latarr = np.array(lat.copy())
    hours = (lat.index - pd.Timestamp('2000-1-1'))/pd.Timedelta('1hour')
    # Find ascending nodes
    pnlat = np.where(latarr>=0,1,-1)
    dpnlat = np.diff(pnlat)
    anode, = np.where(dpnlat == 2)
    # exclude data gap
    dt = hours[anode+1] - hours[anode]
    anode = anode[dt<0.1] #  0.1 hours
    if anode.size < 2:
        return None
    bnode = anode[:-1]+1
    enode = anode[1:]
    for kb, ke in zip(bnode, enode):
        if hours[ke] - hours[kb] > 2:
            continue
        arglat[kb:ke+1] = (360 + latarr[ke] - latarr[kb]) / \
                (hours[ke]-hours[kb]) * \
                (hours[kb:ke+1] - hours[kb]) + latarr[kb]
    return arglat

def great_circle_distance_earth(lat1,lon1,lat2,lon2):
    '''
    Calculate great circle distance: the shortest distance of
    two points in the Earth surface.
      lat1: source latitudes (degree)
      lon1: source longitudes (degree)
      lat2: destination latitudes (degree)
      lon2: destination longitudes (degree)
    '''
    lat1, lon1, lat2, lon2 = (
            np.array(lat1), np.array(lon1), np.array(lat2), np.array(lon2))
    EARTHRADIUS = 6371.009 #  Unit: km
    lat1, lon1 = lat1/180*np.pi, lon1/180*np.pi
    lat2, lon2 = lat2/180*np.pi, lon2/180*np.pi
    dtheta = np.arccos(np.sin(lat1)*np.sin(lat2) +
                       np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    return EARTHRADIUS*dtheta

def set_lat_lt_polar(ax, ns='N', boundinglat=0,
                     latgrids=tuple(np.arange(-90,91,10)),
                     ltgrids=(0,6,12,18)):
    '''
    Set polar coordinates (latitude, local time).
    '''
    fc =1 if ns is 'N' else -1
    rgrids = 90 - fc*np.array(latgrids)
    rgrids = np.sort(rgrids[(rgrids>0) & (rgrids<=90)])
    thetagrids = np.array(ltgrids)/12*180
    ax.set_rgrids(rgrids,['${:d}^\circ$'.format(k) for k in fc*(90-rgrids)])
    ax.set_thetagrids(thetagrids, (thetagrids/180*12).astype(np.int), frac=1.08)
    ax.set_rmax(90-fc*boundinglat)
    ax.set_theta_zero_location('S')
    return

# END
#--------------------------------------------------------------------------------
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import champ_grace as cg
    import goce as gc
    # Test lat2arglat
    #    a = gc.get_goce_data('2009-6-30','2009-12-10 3:0:0')
    #    arglat = lat2arglat(a.lat)
    #    plt.plot(a.index, a.arglat-arglat, 'bo')
    # Test great_circle_distance_earth
    #    lat1, lon1 = 0, 0
    #    lat2 = np.arange(10)
    #    lon2 = np.arange(10)
    #    a = great_circle_distance_earth(lat1, lon1, lat2, lon2)
    # Test set_lat_lt_polar
    ax = plt.subplot(projection='polar')
    set_lat_lt_polar(ax,ns='S',boundinglat=-50)
    plt.show()
