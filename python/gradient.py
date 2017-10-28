#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_gradient() - calculate the gradient of some variable
#-------------------------------------------------------------------------------
import numpy as np
from np import pi

def calc_gradient(var, xyz, coord='spherical'):
    '''
    Calculate the gradient of a variable
    Input: var        = The variable var[nlon, nlat, nalt]
           xyz        = The coordinates
           coord      = The coordinate system
    Output: gradx, grady, gradz
    Note: If coord is 'spherical', xyz should be (lon, lat alt). In a standard
          spherical coordinates (r,theta,phi), theta is southward and 0 theta
          is the north pole, but the given lat is northward and 0 is the equator.
    '''
    if 'sphe' in coord:
        Re = 6371*1000 # Earth radius, unit: m
        lon, lat, alt = xyz

        gradlon = (1.0/((Re+alt)*np.sin(pi/2-lat))) * \
                  np.gradient(var, axis=0) / \
                  np.gradient(lon, axis=0)
        gradlat = (1.0/(Re+alt)) * \
                  np.gradient(var, axis=1) / \
                  np.gradient(lat, axis=1)
        gradr = np.gradient(var, axis=2) / np.gradient(alt, axis=2)

        return gradlon, gradlat, gradr
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
