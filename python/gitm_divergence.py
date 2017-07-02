#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-01-17 08:58, UM
# Comments: Routine to calculate the divergence in a spherical coordinate.
#
# Include: calc_divergence() - Add divergence to GitmBin instance.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt

def calc_divergence(gdata, neuion='neutral', name='divergence'):
    from spacepy.datamodel import dmarray
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
    divergence = (
            1.0/(RR**2)*np.gradient(RR**2*uwind, axis=2) / \
            np.gradient(alt, axis=2) +
            1.0/(RR*np.cos(lat))*(
                np.gradient(nwind*np.cos(lat), axis=1) / \
                np.gradient(lat, axis=1) +
                np.gradient(ewind, axis=0) / np.gradient(lon, axis=0)))
    gdata[name] = dmarray(divergence, attrs={'units':'/s', 'scale':'linear',
                        'name':neuion+' velocity divergence'})
    return


#END
if __name__ == '__main__':
    import gitm_3D_const_alt as g3ca
    from spacepy.datamodel import dmarray
    path = '/home/guod/WD2T/run_imfby/'
    gdata = gitm.GitmBin(path+'run2/data/3DALL_t100323_060000.bin')
    calc_divergence(gdata, 'neu', name='divergence')

    alt = 200
    nlat = 90
    slat = 50
    ax = plt.subplot(polar=True)
    ax, hc = g3ca.contour_single(
            ax, 'divergence', 'polar', gdata, alt=alt,
            nlat=nlat, slat=slat, dlat=10, dlonlt=6,
            lonlt00='S', nzlevels=20, zcolor=None)
    plt.colorbar(hc)
    plt.show()
