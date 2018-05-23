#------------------------------------------------------------------------------
# By Dongjie Guo, 2016-12-06 10:56, UM
# Comments: Routine to make contour or scatter plot of GITM output on a polar
#           or rectangular coordinates as a function of latitude and local
#           time/longitude.
#------------------------------------------------------------------------------

import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import sys


def calculate_centrallon(gdata, plot_type='polar', useLT=True):
    centrallon = 0 # if useLT=False
    if useLT:
        gt = gdata['time']
        ut = (gt.hour + gt.minute/60 + gt.second/3600)
        if 'po' in plot_type:
            centrallon = (0-ut)*15 # polar
        else:
            centrallon = (12-ut)*15 # rectangular
    return centrallon

def contour_data(zkey, gdata, alt=400, ialt=None):
    '''
    Get Latitude, Local Time (or Longitude) and z data for plot
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           alt        = altitude (default 400 km)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt

    Output: Lat, lon, zdata
    ***************************************************************************
    Note: Lon raries in x direction, lat varies in y direction. This is
          different with GITM.
    ***************************************************************************
    '''
    # Find altitude index
    if ialt == None:
        altitude = gdata['Altitude'][0, 0, 2:-2]
        ialt = np.argmin(abs(altitude-alt*1000))+2 # in GITM, unit of alt is m

    # Find the data
    lon0 = np.array(gdata['dLon'][2:-2, 0, 0])
    lat0 = np.array(gdata['dLat'][0, 2:-2, 0])
    zdata0 = np.array(gdata[zkey][2:-2, 2:-2, ialt])
    zdata0, lon0 = add_cyclic_point(zdata0.T, coord=lon0, axis=1)
    lon0, lat0 = np.meshgrid(lon0, lat0)
    return lon0, lat0, zdata0


def vector_data(gdata, species, alt=400, ialt=None):
    '''
    Input: gdata      = gitm bin structure
           species    = 'neutral' or 'ion'
           alt        = altitude (default 400 km)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt

    Output: lat, lt, nwind, ewind
    ***************************************************************************
    Note: Lon raries in x direction, lat varies in y direction. This is
          different with GITM.
    ***************************************************************************
    '''
    # Find altitude index
    if ialt == None:
        altitude = gdata['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m

    # Find the data
    lon0 = np.array(gdata['dLon'][2:-2, 0, 0])
    lat0 = np.array(gdata['dLat'][0, 2:-2, 0])
    if 'neu' in species.lower():
        nwind = np.array(gdata['V!Dn!N (north)'][2:-2, 2:-2, ialt])
        ewind = np.array(gdata['V!Dn!N (east)'][2:-2, 2:-2, ialt])
    elif 'ion' in species.lower():
        nwind = np.array(gdata['V!Di!N (north)'][2:-2, 2:-2, ialt])
        ewind = np.array(gdata['V!Di!N (east)'][2:-2, 2:-2, ialt])
    nwind, lon0 = add_cyclic_point(nwind.T, coord=lon0, axis=1)
    ewind = add_cyclic_point(ewind.T, axis=1)
    lon0, lat0 = np.meshgrid(lon0, lat0)
    return lon0, lat0, ewind, nwind


def convert_vector(lon0, lat0, ewind, nwind, plot_type, projection):
    # This function mainly converts any vector in the (lon, lat) coordinates
    # to other coordinates.
    if 'po' in plot_type:
        if projection is None:
            sys.exit('projection is not given')
        csign = np.ones(lat0.shape)
        csign[lat0<0] = -1
        centrallon = projection.proj4_params['lon_0']
        if projection.proj4_params['lat_0'] == -90:
            centrallon = centrallon-180
        theta = (lon0-centrallon)/180*np.pi
        theta = csign*theta
        xwind = csign*(ewind*np.cos(theta)-nwind*np.sin(theta))
        ywind = csign*(ewind*np.sin(theta)+nwind*np.cos(theta))
        xyz = projection.transform_points(ccrs.PlateCarree(), lon0, lat0)
        x, y = xyz[..., 0], xyz[..., 1]
        return x, y, xwind, ywind
    else:
        xyz = projection.transform_points(ccrs.PlateCarree(), lon0, lat0)
        lon0, lat0 = xyz[..., 0], xyz[..., 1]
        return lon0, lat0, ewind, nwind


def test(g, alt=400, contour=True, zstr='Rho',levels=None, vector=True,
         neuion='neu', scale=500, useLT=True):
    import gitm_create_coordinate as gcc
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure(figsize=[8.38,8.12])
    # Tested parameters
    polar = [True, True, False]
    nrow = [2,2,2]
    ncol = [2,2,1]
    nax = [1,2,2]
    nlat = [90, -30, 90]
    slat = [30, -90, -90]
    dlat = [10, 10, 30]
    for k in range(3):
        ipolar = polar[k]
        pr = 'pol' if ipolar else 'rec'
        centrallon = calculate_centrallon(g, pr, useLT)
        ax, projection = gcc.create_map(
                nrow[k],ncol[k],nax[k], pr, nlat[k], slat[k], centrallon,
                coastlines=False, dlat=dlat[k],
                useLT=useLT, lonticklabel=(1, 1, 1, 1))
        ialt = np.argmin(np.abs(g['Altitude'][0,0,:]/1000-alt))
        # contour
        if contour:
            lon0, lat0, zdata0 = contour_data(zstr, g, ialt=ialt)
            fplat = (lat0[:,0]>=slat[k]) & (lat0[:,0]<=nlat[k])
            lon0, lat0, zdata0 = (k[fplat,:] for k in [lon0, lat0, zdata0])
            hc = ax.contourf(lon0, lat0, zdata0, 21 if levels is None else levels,
                    transform=ccrs.PlateCarree(),cmap='jet', extend='both')
            print(np.max(zdata0), np.min(zdata0))
        # vector
        if vector:
            lon0, lat0, ewind0, nwind0 = vector_data(g,neuion,alt=alt)
            lon0, lat0, ewind0, nwind0 = convert_vector(
                lon0, lat0, ewind0, nwind0, pr, projection)
            hq = ax.quiver(
                lon0,lat0,ewind0,nwind0,scale=scale,scale_units='inches',
                regrid_shape=20)

        # title
        plt.title(zstr+' at '+'%.2f km' % (g['Altitude'][0,0,ialt]/1000),y=1.05)
    plt.show()
    return fig, hc

#END
#------------------------------------------------------------------------------
if __name__ == '__main__':
    pass

