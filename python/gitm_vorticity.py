#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-01-17 08:58, UM
# Comments: Routine to calculate the vorticity in a spherical coordinate.
#
# Include: calc_vorticity() - Add vorticity to GitmBin instance.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt
import gitm_3D_const_alt as g3ca
from spacepy.datamodel import dmarray
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs

def calc_vorticity(gdata, neuion='neutral', name='nvorticity',
                   component='radial'):
    '''
    calculate the vertical vorticity in a spherical coordinate.
    input:
        gdata       : GitmBin structure
        neuion      : for neutral or ion (default: 'neutral')
        name        : pick a name in GitmBin for the vorticity
        component   : 'radial', 'north', 'east'. Direction of the vorticity.
    return:
        add 'nvorticity' or 'ivorticity' to gdata
    '''
    # Calculate vorticity
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
    if 'radial' in component.lower():
        vort = 1.0/(Re+alt)*(
                np.gradient(nwind, axis=0) / \
                (np.cos(lat)*np.gradient(lon, axis=0)) -
                np.gradient(ewind, axis=1) / \
                np.gradient(lat, axis=1) + (ewind*np.tan(lat)))
        vort = vort + 2*(2*np.pi/24/3600)*np.sin(lat) # coriolis effect
    if 'north' in component.lower(): # coriolis effect?
        vort = -1.0/(Re+alt)*(
                1.0/np.cos(lat)*np.gradient(uwind, axis=0) / \
                np.gradient(lon, axis=0) -
                np.gradient((Re+alt)*ewind, axis=2)/np.gradient(alt, axis=2))
        vort = vort + 2*(2*np.pi/24/3600)*np.cos(lat)
    if 'east' in component.lower(): # coriolis effect?
        vort = 1.0/(Re+alt)*(
                np.gradient(-(r+alt)*nwind, axis=2)/np.gradient(alt, axis=2) +
                np.gradient(uwind, axis=1)/np.gradient(lat, axis=1))
    gdata[name] = dmarray(vort, attrs={'units':'m/s^2', 'scale':'linear',
                        'name':neuion+' velocity '+component+' vorticity'})
    return


#END
if __name__ == '__main__':
    fn = '/home/guod/big/raid4/guod/run_imfby/run2c/data/3DALL_t100323_042000.bin'
    gdata = gitm.GitmBin(fn)
    calc_vorticity(gdata, 'neu', name='nvorticity', component='radial')

    alt = 200
    nlat = 90
    slat = 50
    centrallon = g3ca.calculate_centrallon(gdata, 'polar', useLT=True)
    ax, projection = gcc.create_map(
            1, 1, 1, plot_type='polar', nlat=nlat, slat=slat,
            centrallon=centrallon, coastlines=False, dlat=10)
    lon0, lat0, nvort0 = g3ca.contour_data('nvorticity', gdata, alt=alt)
    lon, lat, ewind, nwind = g3ca.vector_data(gdata, 'neu', alt=alt)
    lon, lat, ewind, nwind = g3ca.convert_vector(
            lon, lat, ewind, nwind, 'polar', projection)
    hc = ax.contourf(
            lon0, lat0, nvort0,
            transform=ccrs.PlateCarree(), cmap='seismic')
    ax.quiver(lon, lat, ewind, nwind, regrid_shape=50, scale_units='inches', scale=1000)
    plt.show()
