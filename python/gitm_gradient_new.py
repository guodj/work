#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_gradient() - calculate the gradient
#-------------------------------------------------------------------------------
import gitm
import numpy as np

def calc_gradient_up(g, data):
    """
    data includes ghost cells
    """
    grad = np.ones(g['Altitude'].shape)*np.nan
    Re = 6371*1000 # Earth radius, unit: m
    alt = g['Altitude'][2:-2, 2:-2, 2:-2]
    grad[2:-2, 2:-2, 2:-2] = \
            np.gradient(data[2:-2, 2:-2, 2:-2], axis=2, edge_order=2) \
            / np.gradient(alt, axis=2, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        grad[2:-2, 2, ialt+2] = np.mean(grad[2:-2, 3, ialt+2])
        grad[2:-2, -2, ialt+2] = np.mean(grad[2:-2, -3, ialt+2])
    return grad

def calc_gradient_north(g, data):
    grad = np.ones(g['Altitude'].shape)*np.nan
    Re = 6371*1000 # Earth radius, unit: m
    lat, lon, alt = \
            (g[k][2:-2, 2:-2, 2:-2] \
             for k in ['Latitude', 'Longitude', 'Altitude'])
    grad[2:-2, 2:-2, 2:-2] = \
            (1.0/(Re+alt)) \
            * np.gradient(data[2:-2, 2:-2, 2:-2], axis=1, edge_order=2) \
            / np.gradient(lat, axis=1, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        grad[2:-2, 2, ialt+2] = np.mean(grad[2:-2, 3, ialt+2])
        grad[2:-2, -2, ialt+2] = np.mean(grad[2:-2, -3, ialt+2])
    return grad

def calc_gradient_east(g, data):
    grad = np.ones(g['Altitude'].shape)*np.nan
    Re = 6371*1000 # Earth radius, unit: m
    lat, lon, alt = \
            (g[k][2:-2, 2:-2, 2:-2] \
             for k in ['Latitude', 'Longitude', 'Altitude'])
    grad[2:-2, 2:-2, 2:-2] = \
            (1.0/((Re+alt)*np.cos(lat))) \
            * np.gradient(data[2:-2, 2:-2, 2:-2], axis=0, edge_order=2) \
            / np.gradient(lon, axis=0, edge_order=2)
    for ialt, alt1 in enumerate(g['Altitude'][0, 0, 2:-2]):
        grad[2:-2, 2, ialt+2] = np.mean(grad[2:-2, 3, ialt+2])
        grad[2:-2, -2, ialt+2] = np.mean(grad[2:-2, -3, ialt+2])
    return grad


if __name__ == '__main__':
    pass
    import gitm_create_coordinate as gcc
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import gitm_3D_const_alt as g3ca
    plt.show()
