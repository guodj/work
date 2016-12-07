#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie, 2016-12-06 10:56, UM
# Comments: Routine to make contour or scatter plot of GITM output on a polar or
#           rectangular coordinates.
#
# Include: gitm_3D_lat_lt_single() - a single plot, polar or rectangular, fixed
#                                    altitude, apecified latitudes.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt


def gitm_3D_lat_lt_data(zkey, gdata, alt=400, nlat=90, slat=-90,
                        plot_type='polar', *args, **kwargs):
    '''
    Get Latitude, Local Time and zkey data for scatter or contour plot
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           plot_type  = key to determine plot type (rectangular, polar)
           alt        = altitude (km, default 400 km)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)

    Output: Lat, LT, zdata
    '''
    # Find altitude index
    altitude = gdata['Altitude'][0, 0, :]
    ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
    # Confine latitudes and longitudes
    latitude = gdata['Latitude'][0, :, 0]
    ilat = np.argwhere((latitude >= -np.pi/2) &
                       (latitude <= np.pi/2) &
                       (latitude >= slat*np.pi/180) &
                       (latitude < nlat*np.pi/180))
    ilatmin, ilatmax = ilat.min(), ilat.max()+1
    # max+1 is due to python slice syntax
    longitude = gdata['Longitude'][:, 0, 0]
    ilon = np.argwhere((longitude>0) & (longitude<2*np.pi))
    ilonmin, ilonmax = ilon.min(), ilon.max()+1
    # max+1 is due to python slice syntax
    if 'rec' in plot_type.lower():
        ilatmin, ilatmax = ilatmin-1, ilatmax+1
    # Find the data
    lat0 = gdata['Latitude'][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    lat0 = lat0*180/np.pi
    lt0 = gdata['LT'][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    zdata0 = gdata[zkey][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    # Sort LT
    ilt = np.argmin(lt0[:,0])
    lat0 = np.roll(lat0, -ilt, 0)
    lt0 = np.roll(lt0, -ilt, 0)
    zdata0 = np.roll(zdata0, -ilt, 0)
    # Extend LT (for contourf)
    lat0 = np.insert(lat0, 0, lat0[-1, :], axis=0)
    lat0 = np.insert(lat0[::-1, :], 0, lat0[1, :], axis=0)[::-1, :]
    lt0 = np.insert(lt0, 0, lt0[-1, :], axis=0)
    lt0 = np.insert(lt0[::-1, :], 0, lt0[1, :], axis=0)[::-1, :]
    zdata0 = np.insert(zdata0, 0, zdata0[-1, :], axis=0)
    zdata0 = np.insert(zdata0[::-1, :], 0, zdata0[1, :], axis=0)[::-1, :]
    lt0[0, :] = lt0[0, :]-24
    lt0[-1, :] = lt0[-1, :]+24
    return lat0, lt0, zdata0


def gitm_3D_lat_lt_plot(ax, plot_type, lat0, lt0, zdata0,  title=None,
                        figname=None, draw=True, nlat=90, slat=-90, dlat=30,
                        dlt=6, lt00='S', zmax=None, zmin=None, zcolor=None,
                        data_type="contour", *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: ax         = ax to draw the plot on
           plot_type  = key to determine plot type (rectangular, polar)
           lat0       = 2D latitude data
           lt0        = 2D local time data
           zdata0     = 2D data to be drawn
           title      = plot title
           figname    = file name to save figure as (default is none)
           draw       = draw to screen? (default is True)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 30)
           dlt        = increment of local time ticks (default 6)
           lt00       = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    if (zmin is None) | (zmax is None):
        zmin, zmax = np.min(zdata0), np.max(zdata0)
    # set color map
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic_r"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    if 'pol' in plot_type.lower(): # polar
        if slat*nlat < 0:
            print('Error: For polar plot, nlat and slat should '
                  'have the same sign')
            return
        csign = 1 if nlat > 0 else -1
        # Convert lat and lt to r and theta
        r0 = 90-csign*lat0
        theta0 = lt0/12*np.pi
        hcont = ax.contourf(np.array(theta0), np.array(r0), np.array(zdata0),
                            levels=np.linspace(zmin, zmax, 50),
                            cmap=zcolor, extend='both', *args, **kwargs)
        # Set polar coordinates
        rticks = np.arange(dlat, r0.max()+dlat/2, dlat)
        rlabels = ['{:02.0f}$^\circ$'.format(k) for k in csign*(90-rticks)]
        thetaticks = np.arange(0, 360, dlt*180/12)
        thetalabels = ['{:02.0f}'.format(k) for k in thetaticks*12/180]
        ax.set_rgrids(rticks, rlabels)
        ax.set_thetagrids(thetaticks, thetalabels)
        ax.set_theta_zero_location(lt00)
        return ax, hcont
    if 'rec' in plot_type.lower(): # rectangular
        hcont = ax.contourf(np.array(lt0), np.array(lat0), np.array(zdata0),
                            levels=np.linspace(zmin, zmax, 50),
                            cmap=zcolor, extend='both', *args, **kwargs)
        xticks = np.arange(0, 24+dlt/2, dlt)
        xticklabels = ['{:02.0f}'.format(k) for k in xticks]
        yticks = np.arange(slat, nlat+0.00001, dlat)
        yticklabels = ['{:02.0f}'.format(k) for k in yticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_xlim(0,24)
        ax.set_ylim(slat, nlat)
        return ax, hcont


if __name__ == '__main__':
    import gitm
    import pandas as pd
    fig, ax = plt.subplots(4, 2, figsize=[6, 13],
                           subplot_kw=dict(projection='polar'))
    for k in range(3, 25, 3):
        tt0 = pd.Timestamp('2010-02-26')+pd.Timedelta(k-3, 'h')
        tt1 = pd.Timestamp('2010-02-26')+pd.Timedelta(k, 'h')
        ttstr0 = tt0.strftime('%y%m%d_%H%M%S')
        ttstr1 = tt1.strftime('%y%m%d_%H%M%S')
        a0 = gitm.GitmBin(
                '/home/guod/tmp/3DALL_msis/3DALL_t'+ttstr0+'.bin',
                varlist=['Rho'])
        lat0, lt0, rho0 = gitm_3D_lat_lt_data('Rho', a0, nlat=0, slat=-90)
        a1 = gitm.GitmBin(
                '/home/guod/tmp/3DALL_msis/3DALL_t'+ttstr1+'.bin',
                varlist=['Rho'])
        lat1, lt1, rho1 = gitm_3D_lat_lt_data('Rho', a1, nlat=0, slat=-90)
        plt.sca(ax[int((k/3-1)//2), int(k/3%2-1)])
        axt, hcont = gitm_3D_lat_lt_plot(
                plt.gca(), 'polar', lat0, lt0, 100*(rho1-rho0)/rho0, title=None,
                figname=None, draw=True, nlat=0, slat=-90, dlat=30, dlt=6,
                lt00='S', zmax=15, zmin=-15, zcolor=None, data_type="contour")
        plt.title(tt1.strftime('%H')+'-'+tt0.strftime('%H'))

    [ax[-1, k].set_xlabel('LT') for k in range(2)]
    plt.subplots_adjust(top=0.95, bottom=0.07)
    cax = plt.axes([0.3, 0.025, 0.4, 0.01])
    plt.colorbar(hcont, cax=cax, ticks=np.arange(-15, 16, 5),
                 orientation='horizontal')
    plt.show()
