#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-01-17 08:58, UM
# Comments: Routine to calculate the vertical vorticity in a spherical coordinate.
#
# Include: add_vorticity_r() - Add the radial vorticity (neu or ion) to GitmBin.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt

def add_vorticity_r(gdata, neuion='neutral'):
    from spacepy.datamodel import dmarray
    '''
    calculate the vertical vorticity in a spherical coordinate.
    input:
        gdata       : GitmBin structure
        neuion      : for neutral or ion (default: 'neutral')
    return:
        add 'nvorticity' or 'ivorticity' to gdata
    '''
    # Calculate vorticity
    lat, lon = (gdata[k] for k in ['Latitude', 'Longitude'])
    if 'neu' in neuion.lower():
        nwind = gdata['V!Dn!N (north)']  # unit: m/s
        ewind = gdata['V!Dn!N (east)']  # unit: m/s
    elif 'ion' in neuion.lower():
        nwind = gdata['V!Di!N (north)']  # unit: m/s
        ewind = gdata['V!Di!N (east)']  # unit: m/s
    Re = 6371*1000 # Earth radius, unit: m
    vorticity_r = 1.0/(Re+gdata['Altitude'])*(
            np.gradient(nwind, axis=0) / (np.cos(lat)*np.gradient(lon, axis=0)) -
            np.gradient(ewind, axis=1) / np.gradient(lat, axis=1) + (ewind*np.tan(lat)))
    ac = 'nvorticity' if 'neu' in neuion else 'ivorticity'
    gdata[ac] = dmarray(vorticity_r, attrs={'units':'m/s^2', 'scale':'linear',
                        'name':neuion+' radial vorticity'})
    return


#END
if __name__ == '__main__':
    import gitm_3D_const_alt as g3ca
    from spacepy.datamodel import dmarray
    path = '/home/guod/WD2T/run_imfby/'
    gdata = gitm.GitmBin(path+'run2/data/3DALL_t100323_050000.bin',
                     varlist=['V!Dn!N (north)','V!Dn!N (east)', 'Rho'])
    add_vorticity_r(gdata, 'neu')

    alt = 200
    nlat = -40
    slat = -90
    ax = plt.subplot(polar=True)
    ax, hc = g3ca.contour_single(
        ax, 'nvorticity', 'polar', gdata, alt=alt,
        nlat=nlat, slat=slat, dlat=10, dlonlt=6,
        lonlt00='S', nzlevels=10, zcolor=None)
    ax, hc = g3ca.vector_single(
            ax, gdata, 'neu', 'polar', alt=alt, nlat=nlat, slat=slat, dlat=10,
            dlonlt=6, lonlt00='S',  scale=1000,
            scale_units='inches',  useLT=True)
    plt.show()
