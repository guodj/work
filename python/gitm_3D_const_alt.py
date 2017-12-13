#------------------------------------------------------------------------------
# By Dongjie Guo, 2016-12-06 10:56, UM
# Comments: Routine to make contour or scatter plot of GITM output on a polar
#           or rectangular coordinates as a function of latitude and local
#           time/longitude.
#------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import gitm_create_coordinate as gcc
import sys


def calculate_centrallon(gdata, plot_type='polar', useLT=True):
    centrallon = 0
    if useLT:
        gt = gdata['time']
        ut = (gt.hour + gt.minute/60 + gt.second/3600)
        if 'po' in plot_type:
            centrallon = (0-ut)*15
        else:
            centrallon = (12-ut)*15
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
          different with GITM settings.
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
    return np.array(lon0), np.array(lat0), np.array(zdata0)


def vector_data(gdata, species, alt=400, ialt=None):
    '''
    Obtain data for _vector_plot.
    Input: gdata      = gitm bin structure
           species    = 'neutral' or 'ion'
           alt        = altitude (default 400 km)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt

    Output: lat, lt, nwind, ewind
    ***************************************************************************
    Note: Lon raries in x direction, lat varies in y direction. This is
          different with GITM settings.
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
    return np.array(lon0), np.array(lat0), np.array(ewind), np.array(nwind)


def convert_vector(lon0, lat0, ewind, nwind, plot_type, projection):
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


#END
#------------------------------------------------------------------------------
if __name__ == '__main__':
    # test
    plt.close('all')
    path = '/home/guod/tmp/'
    g = gitm.GitmBin(
            path+'3DMOM_t030322_001000.bin',
            varlist=['NeuPressureGrad (up)', 'GravityForce',
                'VGradV (up)', 'CentriForce (up)', 'CoriolisForce (up)',
                'ViscosityForce (up)','SpheGeomForce (up)'])
    g['forces'] = g['NeuPressureGrad (up)']+g['GravityForce']+g['VGradV (up)']+\
            g['CentriForce (up)']+g['CoriolisForce (up)']+\
            g['SpheGeomForce (up)']+g['ViscosityForce (up)']

    # Tested parameters
    polar = True
    useLT = True
    nlat, slat, dlat, dlon = 90, 0, 10, 90
    #-----------------------
    pr = 'pol' if polar else 'rec'
    centrallon = calculate_centrallon(g, pr, useLT)
    ax, projection = gcc.create_map(
            1,1,1, pr, nlat, slat, centrallon, coastlines=False, dlat=dlat,
            dlon=dlon, useLT=useLT, lonticklabel=(1, 1, 1, 1), aspect=1.5)
    lon0, lat0, zdata0 = contour_data('forces', g, alt=400)
    hc = ax.contourf(lon0, lat0, zdata0, transform=ccrs.PlateCarree(),cmap='jet')
    plt.colorbar(hc)
    # lon0, lat0, ewind, nwind = vector_data(g, 'neutral', alt=400)
    # lon0, lat0, ewind, nwind = convert_vector(
    #         lon0, lat0, ewind, nwind, plot_type=pr, projection=projection)
    # ax.quiver(lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
    #           regrid_shape=30)
    plt.show()
