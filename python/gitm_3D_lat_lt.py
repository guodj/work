#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie, 2016-12-06 10:56, UM
# Comments: Routine to make contour or scatter plot of GITM output on a polar or
#           rectangular coordinates as a function of (magnetic) latitude and
#           (magnetic) local time.
#
# Include: gitm_3D_lat_lt_data() - obtain 2D latitude, local time and z data for
#                  a contour or scatter plot. Fixed altitude, specified
#                  latitudes.
#          gitm_3D_lat_lt_plot() - contour or scatter plot on a polar or
#                  rectangular coordinates. Fixed altitude, specified latitudes.
#          gitm_3D_lat_lt_single() - contour or scatter plot on a polar or
#                  rectangular coordinates. Fixed altitude, specified latitudes.
#                  For 1 Gitm output file.
#          gitm_3D_lat_lt_diff() - contour or scatter plot on a polar or
#                  rectangular coordinates. Fixed altitude, specified latitudes.
#                  For difference of 2 Gitm output files.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt


def gitm_3D_lat_lt_data(zkey, gdata, alt=400, nlat=90, slat=-90):
    '''
    Get Latitude, Local Time and z data for scatter or contour plot
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           alt        = altitude (default 400 km)
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
    ilatmin, ilatmax = ilat.min(), ilat.max()
    longitude = gdata['Longitude'][:, 0, 0]
    ilon = np.argwhere((longitude>0) & (longitude<2*np.pi))
    ilonmin, ilonmax = ilon.min(), ilon.max()
    # Find the data
    lat0 = gdata['Latitude'][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    lat0 = lat0*180/np.pi
    lt0 = gdata['LT'][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    zdata0 = gdata[zkey][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    # Sort LT
    ilt = np.argmin(lt0[:,0]) # LT range [0, 24]
    lat0 = np.roll(lat0, -ilt, 0)
    lt0 = np.roll(lt0, -ilt, 0)
    zdata0 = np.roll(zdata0, -ilt, 0)
    # Extend latitude
    if nlat == 90:
        lat0 = np.insert(lat0[:, ::-1], 0, 90, axis=1)[:, ::-1]
        lt0 = np.insert(lt0[:, ::-1], 0, lt0[:, -1], axis=1)[:, ::-1]
        zdata0 = np.insert(
                zdata0[:, ::-1], 0, np.mean(zdata0[:, -1]), axis=1)[:, ::-1]
    if slat == -90:
        lat0 = np.insert(lat0, 0, -90, axis=1)
        lt0 = np.insert(lt0, 0, lt0[:,0], axis=1)
        zdata0 = np.insert(zdata0, 0, np.mean(zdata0[:, 0]), axis=1)
    # Extend LT (for contour)
    lat0 = np.insert(lat0, 0, lat0[-1, :], axis=0)
    lat0 = np.insert(lat0[::-1, :], 0, lat0[1, :], axis=0)[::-1, :]
    lt0 = np.insert(lt0, 0, lt0[-1, :], axis=0)
    lt0 = np.insert(lt0[::-1, :], 0, lt0[1, :], axis=0)[::-1, :]
    zdata0 = np.insert(zdata0, 0, zdata0[-1, :], axis=0)
    zdata0 = np.insert(zdata0[::-1, :], 0, zdata0[1, :], axis=0)[::-1, :]
    lt0[0, :] = lt0[0, :]-24
    lt0[-1, :] = lt0[-1, :]+24
    return lat0, lt0, zdata0


def gitm_3D_mlat_mlt_data(zkey, gdata, alt=400, nlat=90, slat=-90):
    '''
    Get magnetic Latitude, Local Time and z data for scatter or contour plot
    Input: zkey        = key for z variable (ie 'Vertical TEC')
           gdata       = gitm bin structure
           alt         = altitude (default 400 km)
           nlat        = northern magnetic latitude limit
                         (degrees North, default 90)
           slat        = southern magneticlatitude limit
                         (degrees North, defalut -90)

    Output: mLat, mLT, zdata
    '''
    from apexpy import Apex
    from scipy.spatial import cKDTree
    # Find altitude index
    altitude = gdata['Altitude'][0, 0, :]
    ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
    # Confine latitudes and longitudes
    latitude = gdata['Latitude'][0, :, 0]
    ilat = np.argwhere((latitude >= -np.pi/2) &
                       (latitude <= np.pi/2))
    ilatmin, ilatmax = ilat.min(), ilat.max()
    longitude = gdata['Longitude'][:, 0, 0]
    ilon = np.argwhere((longitude>0) & (longitude<2*np.pi))
    ilonmin, ilonmax = ilon.min(), ilon.max()
    lat0 = gdata['Latitude'][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    lat0 = lat0*180/np.pi
    lon0 = gdata['Longitude'][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    lon0 = lon0*180/np.pi
    lt0 = gdata['LT'][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    zdata0 = gdata[zkey][ilonmin:ilonmax+1, ilatmin:ilatmax+1, ialt]
    # Convert lat, lt to mlat, mlt
    yearf = gdata['time'].year+gdata['time'].month/12
    apx = Apex(date=yearf)
    mlat0, mlt0 = apx.convert(
            lat0, lon0, source='geo', dest='mlt',
            height=gdata['Altitude'][0, 0, ialt]/1000,
            datetime=gdata['time'])
    mlat0, mlt0 = mlat0.reshape(-1, 1), mlt0.reshape(-1, 1)
    zdata0 = zdata0.reshape(-1, 1)
    # Convert spherical coordinate to cartesion coordinate
    x0 = np.cos(mlat0*np.pi/180)*np.sin(mlt0*np.pi/12)
    y0 = np.cos(mlat0*np.pi/180)*np.cos(mlt0*np.pi/12)
    z0 = np.sin(mlat0*np.pi/180)
    # Create cKDtree
    tree = cKDTree(np.concatenate([x0, y0, z0], axis=1))
    # Create magnetic coordinate grids
    mlat1 = np.meshgrid(np.arange(slat, nlat+1, 1))
    mlt1 = np.meshgrid(np.arange(0, 24.000001, 3/60))
    mlat1, mlt1 = np.meshgrid(mlat1, mlt1)
    # Convert spherical coordinate to cartesion coordinate
    x1 = np.cos(mlat1*np.pi/180)*np.sin(mlt1*np.pi/12)
    y1 = np.cos(mlat1*np.pi/180)*np.cos(mlt1*np.pi/12)
    z1 = np.sin(mlat1*np.pi/180)
    # Get results
    x1, y1, z1 = x1.reshape(-1, 1), y1.reshape(-1, 1), z1.reshape(-1, 1)
    d, izdata1 = tree.query(np.concatenate([x1, y1, z1], axis=1), k=1)
    #w = 1.0/d**2
    #zdata1 = np.sum(
    #        w*(zdata0[izdata1].reshape(w.shape)), axis=1)/np.sum(w, axis=1)
    #zdata1 = zdata1.reshape(mlat1.shape)
    zdata1 = zdata0[izdata1].reshape(mlat1.shape)
    return  mlat1, mlt1, zdata1


def gitm_3D_lat_lt_plot(ax, plot_type, lat0, lt0, zdata0,
                        nlat=90, slat=-90, dlat=10, dlt=6, lt00='S',
                        zmax=None, zmin=None, zcolor=None,
                        data_type="contour", *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: ax         = ax to draw the plot on
           plot_type  = key to determine plot type (rectangular, polar)
           lat0       = 2D latitude data
           lt0        = 2D local time data
           zdata0     = 2D data to be drawn
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
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
        if 'cont' in data_type.lower(): # contour
            hcont = ax.contourf(np.array(theta0), np.array(r0),
                                np.array(zdata0),
                                levels=np.linspace(zmin, zmax, 21),
                                cmap=zcolor, extend='both', *args, **kwargs)
        if 'sca' in data_type.lower(): # scatter
            hcont = ax.scatter(np.array(theta0), np.array(r0),
                               c=np.array(zdata0), vmin=zmin, vmax=zmax,
                               cmap=zcolor, lw=0, *args, **kwargs)
        # Set polar coordinates
        rticks = np.arange(dlat, r0.max()+dlat/2, dlat)
        rlabels = ['{:02.0f}$^\circ$'.format(k) for k in csign*(90-rticks)]
        thetaticks = np.arange(0, 360, dlt*180/12)
        thetalabels = ['{:02.0f}'.format(k) for k in thetaticks*12/180]
        ax.set_rgrids(rticks, rlabels)
        ax.set_thetagrids(thetaticks, thetalabels)
        ax.set_theta_zero_location(lt00)
        ax.set_rlim(0, min(90-slat, 90+nlat))
        return ax, hcont
    if 'rec' in plot_type.lower(): # rectangular
        if 'cont' in data_type.lower(): # contour
            hcont = ax.contourf(np.array(lt0), np.array(lat0), np.array(zdata0),
                                levels=np.linspace(zmin, zmax, 21),
                                cmap=zcolor, extend='both', *args, **kwargs)
        if 'sca' in data_type.lower(): # scatter
            hcont = ax.scatter(np.array(lt0), np.array(lat0),
                               c=np.array(zdata0), vmin=zmin, vmax=zmax,
                               cmap=zcolor, lw=0, *args, **kwargs)
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


def gitm_3D_lat_lt_single(ax, zkey, plot_type, gdata, alt=400, title=False,
                          nlat=90, slat=-90, dlat=10, dlt=6, mag=False,
                          lt00='S', zmax=None, zmin=None, zcolor=None,
                          data_type="contour", *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range at a fixed altitude for one 3D file.
    Input: ax         = ax to draw the plot on
           zkey       = key for z variable (ie 'Vertical TEC')
           plot_type  = key to determine plot type (rectangular, polar)
           gdata      = gitm bin structure
           alt        = altitude (km, default 400 km)
           title      = whether to use default title (default False)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
           dlt        = increment of local time ticks (default 6)
           mag        = whether use magnetic coordinates (default False)
           lt00       = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    from apexpy import Apex
    if not mag:
        lat0, lt0, zdata0 = gitm_3D_lat_lt_data(
                zkey, gdata, alt=alt, nlat=nlat, slat=slat)
    if mag:
        lat0, lt0, zdata0 = gitm_3D_mlat_mlt_data(
                zkey, gdata, alt=alt, nlat=nlat, slat=slat)
    ax, hc = gitm_3D_lat_lt_plot(
             ax, plot_type, lat0, lt0, zdata0,
             nlat=nlat, slat=slat, dlat=dlat,
             dlt=dlt, lt00=lt00, zmax=zmax, zmin=zmin, zcolor=zcolor,
             data_type=data_type, *args, **kwargs)
    if title:
        # Find altitude index
        altitude = gdata['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(gdata['time'].strftime('%y-%m-%d %H:%M:%S')+
                     ' @ '+'%5.1f'%(altalt/1000)+' km')
    return ax, hc


def gitm_3D_lat_lt_diff(
        ax, zkey, plot_type, gdata1, gdata2, alt=400, diff_type='relative',
        title=False,  nlat=90, slat=-90, dlat=10, dlt=6, mag=False,
        lt00='S', zmax=None, zmin=None, zcolor=None, data_type="contour",
        *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range at a fixed altitude for one 3D file.
    Input: ax         = ax to draw the plot on
           zkey       = key for z variable (ie 'Vertical TEC')
           plot_type  = key to determine plot type (rectangular, polar)
           gdata1     = gitm bin structure 1
           gdata2     = gitm bin structure 2
           alt        = altitude (km, default 400 km)
           diff_type  = relative or absolute density difference
                        (default relative)
           title      = whether to use default title (default False)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
           dlt        = increment of local time ticks (default 6)
           mag        = whether use magnetic coordinates (default False)
           lt00       = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    if not mag:
        lat1, lt1, zdata1 = gitm_3D_lat_lt_data(
                zkey, gdata1, alt=alt, nlat=nlat, slat=slat)
        lat2, lt2, zdata2 = gitm_3D_lat_lt_data(
                zkey, gdata2, alt=alt, nlat=nlat, slat=slat)
    if mag:
        lat1, lt1, zdata1 = gitm_3D_mlat_mlt_data(
                zkey, gdata1, alt=alt, nlat=nlat, slat=slat)
        lat2, lt2, zdata2 = gitm_3D_mlat_mlt_data(
                zkey, gdata2, alt=alt, nlat=nlat, slat=slat)
    if 'rel' in diff_type.lower():
        zdata = 100*(zdata2-zdata1)/zdata1
    if 'abs' in diff_type.lower():
        zdata = zdata2-zdata1
    ax, hc =gitm_3D_lat_lt_plot(
            ax, plot_type, lat1, lt1, zdata,
            nlat=nlat, slat=slat, dlat=dlat,
            dlt=dlt, lt00=lt00, zmax=zmax, zmin=zmin, zcolor=zcolor,
            data_type=data_type, *args, **kwargs)
    if title:
        # Find altitude index
        altitude = gdata1['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(gdata2['time'].strftime('%m-%d %H:%M')+' - '+
                     gdata1['time'].strftime('%m-%d %H:%M')+' @ '+
                     '%5.1f'%(altalt/1000)+' km')
    return ax, hc


if __name__ == '__main__':
    # Test gitm_3D_lat_lt_single and gitm_3D_lat_lt_diff
    import gitm
    import pandas as pd
    path = '/home/guod/WD2T/run_imfby/'
    g1 = gitm.GitmBin(
            path+'run1/data/3DALL_t100323_010000.bin', varlist=['Rho'])
    g2 = gitm.GitmBin(
            path+'run2/data/3DALL_t100323_010000.bin', varlist=['Rho'])
    ax = plt.subplot(polar=True)
    ax, hc = gitm_3D_lat_lt_diff(
            ax, 'Rho', 'pol', g1, g2, alt=400, diff_type='relative',
            title=True, nlat=90, slat=40, dlat=10,
            dlt=6, lt00='S', zmax=20, zmin=-20, zcolor=None, mag=False,
            data_type="contour")
    plt.colorbar(hc)
    plt.show()
