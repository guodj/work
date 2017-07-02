#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_gradient() - calculate the gradient of some variable in gitm
#-------------------------------------------------------------------------------
import gitm
import numpy as np

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
    from spacepy.datamodel import dmarray
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
    g = gitm.GitmBin(
            '/home/guod/WD2T/run_imfby/run1/data/3DALL_t100323_012500.bin')
    calc_gradient(g, 'Rho', 'gradient_r', component='radial')
