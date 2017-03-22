#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-01-17 08:58, UM
# Comments: Routine to calculate the vorticity in a spherical coordinate.
#
# Include: calc_vorticity() - Add vorticity to GitmBin instance.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt

def calc_vorticity(gdata, neuion='neutral', name='nvorticity', component='radial'):
    from spacepy.datamodel import dmarray
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
                np.gradient(nwind, axis=0) / (np.cos(lat)*np.gradient(lon, axis=0)) -
                np.gradient(ewind, axis=1) / np.gradient(lat, axis=1) + (ewind*np.tan(lat)))
    if 'north' in component.lower():
        vort = -1.0/(Re+alt)*(
                1.0/np.cos(lat)*np.gradient(uwind, axis=0)/np.gradient(lon, axis=0) -
                np.gradient((Re+alt)*ewind, axis=2)/np.gradient(alt, axis=2))
    if 'east' in component.lower():
        vort = 1.0/(Re+alt)*(
                np.gradient(-(r+alt)*nwind, axis=2)/np.gradient(alt, axis=2) +
                np.gradient(uwind, axis=1)/np.gradient(lat, axis=1))
    gdata[name] = dmarray(vort, attrs={'units':'m/s^2', 'scale':'linear',
                        'name':neuion+' velocity '+component+' vorticity'})
    return


#END
if __name__ == '__main__':
    import gitm_3D_const_alt as g3ca
    from spacepy.datamodel import dmarray
    path = '/home/guod/WD2T/run_imfby/'
    gdata = gitm.GitmBin(path+'run2/data/3DALL_t100323_000000.bin')
    calc_vorticity(gdata, 'neu', name='nvorticity')

    alt = 200
    nlat = 90
    slat = 50
    ax = plt.subplot(polar=True)
    ax, hc = g3ca.contour_single(
            ax, 'nvorticity', 'polar', gdata, alt=alt,
            nlat=nlat, slat=slat, dlat=10, dlonlt=6,
            lonlt00='S', nzlevels=20, zcolor=None)
    ax, hc = g3ca.vector_single(
            ax, gdata, 'neu', 'polar', alt=alt, nlat=nlat, slat=slat, dlat=10,
            dlonlt=6, lonlt00='S',  scale=1000,
            scale_units='inches',  useLT=True)
    plt.show()
