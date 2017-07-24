#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_gradient() - calculate the gradient of some variable in gitm
#-------------------------------------------------------------------------------
import gitm
import numpy as np
from spacepy.datamodel import dmarray
import gitm_create_coordinate as gcc
import gitm_3D_const_alt as g3ca
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def calc_gradient(g, var, name, component='radial'):
    '''
    Calculate the gradient of a variable in gitm.
    Input: g          = gitm data
           var        = variable used to calculate gradient
           component  = 'radial', 'north' or 'east'.
                        On which direction the gradient is
                        calculated (default 'radial')
           name       = pick a name for the gradient of `var`.
    Output: g with `name` added
    '''
    Re = 6371*1000 # Earth radius, unit: m
    lat, lon, alt, zdata = (g[k] for k in
                            ['Latitude', 'Longitude', 'Altitude', var])
    if 'radial' in component.lower():
        grad = np.gradient(zdata, axis=2) / np.gradient(alt, axis=2)
    if 'north' in component.lower():
        grad = (1.0/(Re+alt))*np.gradient(zdata, axis=1) / \
                np.gradient(lat, axis=1)
    if 'east' in component.lower():
        grad = ((1.0/((Re+alt)*np.sin(np.pi/2-lat))) *
                np.gradient(zdata, axis=0) / np.gradient(lon, axis=0))
    g[name] = dmarray(
            grad,
            attrs={'units':'None', 'scale':'linear',
                   'name':(g[var].attrs['name']+' '+component+
                           ' component gradient')})
    return g


if __name__ == '__main__':
    plt.close('all')
    g = gitm.GitmBin(
            '/home/guod/big/raid4/guod/run_imfby/run1c/data/'
            '3DALL_t100323_012500.bin')
    calc_gradient(g, 'Rho', 'gradn', component='north')
    calc_gradient(g, 'Rho', 'grade', component='east')
    centrallon = g3ca.calculate_centrallon(g, 'polar', useLT=True)
    ax, projection = gcc.create_map(
            1, 1, 1, plot_type='polar', nlat=90, slat=50,
            centrallon=centrallon, coastlines=False, dlat=10)
    lon, lat, rho = g3ca.contour_data('Rho', g)
    ax.contourf(lon, lat, rho, transform=ccrs.PlateCarree(),
                levels=np.linspace(1e-12, 7.5e-12, 21), cmap='viridis')
    lon, lat, gradn = g3ca.contour_data('gradn', g)
    lon, lat, grade = g3ca.contour_data('grade', g)
    lon, lat, ewind, nwind = g3ca.convert_vector(
            lon, lat, grade, gradn, plot_type='polar', projection=projection)
    ax.quiver(lon, lat, ewind*1e21, nwind*1e21, regrid_shape=50)
    plt.show()
