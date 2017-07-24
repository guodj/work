#--------------------------------------------------------------------------------
# Functions that may be used many times
# By Dongjie, start from Mon Oct 10 09:37:12 EDT 2016
#--------------------------------------------------------------------------------

# Global imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import pi

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

def set_polar(ax, ns='N', boundinglat=0, dlat=10, dlon=90, useLT=True):
    '''
    Set polar coordinate
    Input:
        ax          : on which the polar coordinate is created
        ns          : 'N' or 'S', northern or southern hemispheres
        boundinglat : latitude at the boundary
        dlat        : latitude grid gap, unit is degree, positive even for SH
        dlon        : longitude grid gap, unit is degree even useLT is True, positive
        useLT       : if True, longitude label is LT
    Return:
        ax
    '''
    if ns=='N':
        rmax = 90-boundinglat
        rtick = np.arange(dlat, rmax+1e-10, dlat)
        rtickl = ['{:.0f}\u00b0'.format(90-k) for k in rtick]
        thetatick = np.arange(0, 360, dlon)
        thetatickl = ['{:.0f}\u00b0'.format(k) for k in thetatick]
        if useLT:
            thetatickl = ['{:02.0f}'.format(k*12/180) for k in thetatick]
    if ns=='S':
        rmax = 90+boundinglat
        rtick = np.arange(dlat, rmax+1e-10, dlat)
        rtickl = ['{:.0f}\u00b0'.format(k-90) for k in rtick]
        thetatick = np.arange(0, 360, dlon)
        thetatickl = ['{:.0f}\u00b0'.format((360-k)%360) for k in thetatick]
        if useLT:
            thetatickl = ['{:02.0f}'.format(((360-k)%360)*12/180) for k in thetatick]
    ax.set_rlim(0, rmax)
    ax.set_theta_zero_location('S')
    ax.set_rgrids(rtick, rtickl)
    ax.set_thetagrids(thetatick, thetatickl)
    return ax


def declination_of_sun(doy):
    '''
    note: doy is from 0 (January 1 == 0) and can include decimals
    declination range [-23.44, 23.44]
    '''
    return -np.arcsin(
            0.39779 *
            np.cos(0.98565/180*pi*(doy+10) +
                   1.914/180*pi*np.sin(0.98565/180*pi*(doy-2))))

def solar_zenith_angle(doy, lat, lt):
    return 180/pi*np.arccos(
            np.sin(lat/180*pi)*np.sin(declination_of_sun(doy)) +
            np.cos(lat/180*pi)*np.cos(declination_of_sun(doy)) *
            np.cos((lt-12)/12*pi))


def subplots_create_cbar_axis(ax, location='right', pad=0.02, size=0.05):
    axlc = ax.get_position()
    if location=='right':
        x0 = axlc.x0+axlc.width+pad*axlc.width
        y0 = axlc.y0
        width = size*axlc.width
        height = axlc.height
    elif location=='top':
        x0 = axlc.x0
        y0 = axlc.y0+axlc.height+pad*axlc.height
        width = axlc.width
        height = size*axlc.height
    elif location=='bottom':
        x0 = axlc.x0
        y0 = axlc.y0-pad*axlc.height
        width = axlc.width
        height = size*axlc.height
    return plt.axes([x0, y0, width, height])


def subplots_xylabel(ax, nrow, ncol, xlabel, ylabel, mode='L', **kwargs):
    if type(xlabel)==str:
        xlabel = ncol*[xlabel]
    if type(ylabel)==str:
        ylabel = nrow*[ylabel]
    if mode=='L':
        if nrow == 1:
            ax[0].set_ylabel(ylabel[0], **kwargs)
            for k00 in range(ncol):
                ax[k00].set_xlabel(xlabel[k00], **kwargs)
        elif ncol == 1:
            ax[-1].set_xlabel(xlabel[0], **kwargs)
            for k00 in range(nrow):
                ax[k00].set_ylabel(ylabel[k00], **kwargs)
        else:
            for k00 in range(nrow):
                ax[k00, 0].set_ylabel(ylabel[k00], **kwargs)
            for k00 in range(ncol):
                ax[-1, k00].set_xlabel(xlabel[k00], **kwargs)
    if mode=='1':
        if nrow == 1:
            axt = ax[0]
        elif ncol == 1:
            axt = ax[-1]
        else:
            axt = ax[-1, 0]
        axt.set_ylabel(ylabel[0], **kwargs)
        axt.set_xlabel(xlabel[0], **kwargs)
    return


def subplots_create_abcde(ax, direction='row', x=0.01, y=1.05):
    abcde = list('abcdefghijklmnopqrstuvwxyz')
    if ax.size == ax.shape[0]:
        for k11, k1 in enumerate(ax):
            k1.text(x, y, '( '+abcde[k11]+' )', transform=k1.transAxes)
    else:
        for k1 in range(ax.shape[0]):
            for k2 in range(ax.shape[1]):
                if 'row' in direction:
                    idx=k2+k1*ax.shape[1]
                else:
                    idx=k1+k2*ax.shape[0]
                ax[k1, k2].text(x, y, '( '+abcde[idx]+' )',
                                transform=ax[k1, k2].transAxes)
    return

# END
#-------------------------------------------------------------------------------
if __name__ == '__main__':
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
    # Test set_polar
    #    ax = plt.subplot(polar=True)
    #    ax = set_polar(ax, ns='S', boundinglat=0, dlat=10, dlon=90, useLT=False)
    #    plt.show()

    lat=np.arange(-90, 90)
    lt = np.arange(0, 24)
    doy = np.arange(1, 365)
    plt.plot(lat, solar_zenith_angle(90, lat, 12))
    #plt.plot(doy, declination_of_sun(doy)/np.pi*180)
    plt.show()
