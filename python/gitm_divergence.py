#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-01-17 08:58, UM
# Comments: Routine to calculate the divergence in a spherical coordinate.
#
# Include: calc_divergence() - Add divergence to GitmBin instance.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt
import gitm_3D_const_alt as g3ca
from spacepy.datamodel import dmarray
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs

def calc_divergence(gdata, neuion='neutral', name='divergence'):
    '''
    calculate the divergence in a spherical coordinate.
    input:
        gdata       : GitmBin structure
        neuion      : for neutral or ion (default: 'neutral')
        name        : pick a name in GitmBin for the divergence
    return:
        add neutral or ion divergence to gdata
    Note:
        polar angle = 90-latitude
        azimuthal = longitude
        in spherical delta theta = - delta latitude
    '''
    # Calculate divergence
    lat, lon, alt = (gdata[k] for k in ['Latitude', 'Longitude', 'Altitude'])
    if 'neu' in neuion.lower():
        nwind = gdata['V!Dn!N (north)']  # unit: m/s
        ewind = gdata['V!Dn!N (east)']  # unit: m/s
        uwind = gdata['V!Dn!N (up)']
    elif 'ion' in neuion.lower():
        nwind = gdata['V!Di!N (north)']  # unit: m/s
        ewind = gdata['V!Di!N (east)']  # unit: m/s
        uwind = gdata['V!Di!N (up)']
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt
    ewind = ewind + ((2*np.pi)/(24*3600))*RR*np.cos(lat)
    divergence = (
            1.0/(RR**2)*np.gradient(RR**2*uwind, axis=2) / \
            np.gradient(alt, axis=2) +
            1.0/(RR*np.cos(lat))*(
                np.gradient(nwind*np.cos(lat), axis=1) / \
                np.gradient(lat, axis=1) +
                np.gradient(ewind, axis=0) / np.gradient(lon, axis=0)))
    gdata[name] = dmarray(
            divergence, attrs={'units':'/s', 'scale':'linear',
            'name':neuion+' velocity divergence'})
    return


#END
if __name__ == '__main__':
    plt.close('all')
    path = '/home/guod/big/raid4/guod/run_imfby/run2c/data/'
    gdata = gitm.GitmBin(path+'3DALL_t100323_060000.bin')
    calc_divergence(gdata, 'neu', name='divergence')

    alt = 400
    nlat = 90
    slat = 50
    centrallon = g3ca.calculate_centrallon(gdata, 'polar', useLT=True)
    ax, projection = gcc.create_map(
            1, 1, 1, plot_type='polar', nlat=nlat, slat=slat,
            centrallon=centrallon, coastlines=False, dlat=10)
    lon0, lat0, div0 = g3ca.contour_data('divergence', gdata, alt=alt)
    lon, lat, ewind, nwind = g3ca.vector_data(gdata, 'neu', alt=alt)
    lon, lat, ewind, nwind = g3ca.convert_vector(
            lon, lat, ewind, nwind, 'polar', projection)
    hc = ax.contourf(
            lon0, lat0, div0, levels=np.linspace(-2e-4, 2e-4, 21),
            transform=ccrs.PlateCarree(), cmap='seismic')
    ax.quiver(lon, lat, ewind, nwind, regrid_shape=20)
    plt.show()
