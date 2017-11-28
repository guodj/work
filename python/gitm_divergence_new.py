#-------------------------------------------------------------------------------
# By Dongjie Guo, 2017-11-27 11:17, USTC
# Comments: Routine to calculate the divergence in a spherical coordinate.
#-------------------------------------------------------------------------------

import gitm
import numpy as np

def calc_divergence_up(g, datar):
    divr = np.ones(g['Altitude'].shape)*np.nan
    lat, lon, alt = (g[k][2:-2, 2:-2, 2:-2]
                     for k in ['Latitude', 'Longitude', 'Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt
    RR2 = RR**2
    divr[2:-2, 2:-2, 2:-2] = \
            1.0/RR2 \
            * np.gradient(RR2*datar[2:-2, 2:-2, 2:-2], axis=2, edge_order=2) \
            / np.gradient(alt, axis=2, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        divr[2:-2, 2, ialt+2] = np.mean(divr[2:-2, 3, ialt+2])
        divr[2:-2, -2, ialt+2] = np.mean(divr[2:-2, -3, ialt+2])
    return divr


def calc_divergence_north(g, datan):
    divn = np.ones(g['Altitude'].shape)*np.nan
    lat, lon, alt = (g[k][2:-2, 2:-2, 2:-2]
                     for k in ['Latitude', 'Longitude', 'Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt
    divn[2:-2, 2:-2, 2:-2] = \
            1.0/(RR*np.cos(lat)) \
            * np.gradient(datan[2:-2, 2:-2, 2:-2]*np.cos(lat), axis=1, edge_order=2) \
            / np.gradient(lat, axis=1, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        divn[2:-2, 2, ialt+2] = np.mean(divn[2:-2, 3, ialt+2])
        divn[2:-2, -2, ialt+2] = np.mean(divn[2:-2, -3, ialt+2])
    return divn


def calc_divergence_east(g, datae):
    dive = np.ones(g['Altitude'].shape)*np.nan
    lat, lon, alt = (g[k][2:-2, 2:-2, 2:-2]
                     for k in ['Latitude', 'Longitude', 'Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt
    dive[2:-2, 2:-2, 2:-2] = \
            1.0/(RR*np.cos(lat)) \
            * np.gradient(datae[2:-2, 2:-2, 2:-2], axis=0, edge_order=2) \
            / np.gradient(lon, axis=0, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        dive[2:-2, 2, ialt+2] = np.mean(dive[2:-2, 3, ialt+2])
        dive[2:-2, -2, ialt+2] = np.mean(dive[2:-2, -3, ialt+2])
    return dive


#END
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import gitm_3D_const_alt as g3ca
    from spacepy.datamodel import dmarray
    import gitm_create_coordinate as gcc
    import cartopy.crs as ccrs
    pass
