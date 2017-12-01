import numpy as np

def limiter_mc(dup, ddown, betalimiter=2):
    fp1 = dup>0
    fp2 = ddown>0
    o1 = np.where(dup<ddown, betalimiter*dup,betalimiter*ddown)
    o1 = np.where(o1<(dup+ddown)*0.5, o1, (dup+ddown)*0.5)

    o2 = np.where(dup>ddown, betalimiter*dup,betalimiter*ddown)
    o2 = np.where(o2>(dup+ddown)*0.5, o2, (dup+ddown)*0.5)

    lmc = np.zeros(dup.shape)
    lmc[(fp1) & (fp2)] = o1[(fp1) & (fp2)]
    lmc[(~fp1) & (~fp2)] = o2[(~fp1) & (~fp2)]

    return lmc

def calc_rusanov_lats(lats, alts, var):
    """
    lats: three dimensional with ghost cells (lon, lat,alt)
    alts: three dimensional with ghost cells (lon, lat,alt)
    var: three dimensional with ghost cells (lon, lat,alt)
    """
    nlons = var.shape[0]
    nlats = var.shape[1]
    nalts = var.shape[2]
    dlatdist_gb = np.ones(var.shape)
    dlatdist_fb = np.ones(var.shape)
    dvarlimited = np.ones([nlons, nlats-2, nalts])
    varleft = np.ones([nlons, nlats-3,nalts])
    varright = np.ones([nlons, nlats-3,nalts])
    gradvar = np.ones(var.shape)*np.nan

    Re = 6371*1000
    rdist = Re+alts
    difflat = np.diff(lats, axis=1)


    dlatdist_gb[1:-1, 1:-1,:] = \
            0.5 * (difflat[1:-1,1:,:] + difflat[1:-1,:-1,:])\
            * (rdist[1:-1, 1:-1,:])
    dlatdist_gb[:, 0,:] = dlatdist_gb[:,1,:]
    dlatdist_gb[:,-1,:] = dlatdist_gb[:,-2,:]
    dlatdist_gb[0, :,:] = dlatdist_gb[1,:,:]
    dlatdist_gb[-1,:,:] = dlatdist_gb[-2,:,:]

    dlatdist_fb[1:-1, 1:-1,:] = \
            (difflat[1:-1, :-1,:]) \
            * 0.5 * (rdist[1:-1, 1:-1,:]+ rdist[1:-1, 0:-2,:])
    dlatdist_fb[:, 0,:] = dlatdist_fb[:,1,:]
    dlatdist_fb[:,-1,:] = dlatdist_fb[:,-2,:]
    dlatdist_fb[0, :,:] = dlatdist_fb[1,:,:]
    dlatdist_fb[-1,:,:] = dlatdist_fb[-2,:,:]

    invdlatdist_gb = 1.0/dlatdist_gb
    invdlatdist_fb = 1.0/dlatdist_fb

    factor1 = 0.6250000
    factor2 = 0.0416667

    i = 0
    h = invdlatdist_fb[:,i+2,:]*2
    dvarup = h*(factor1*(var[:,i+2,:]-var[:,i+1,:])
              - factor2*(var[:,i+3,:]-var[:,i,:]))
    dvardown = (var[:, i+1,:] -var[:, i,:])*invdlatdist_fb[:, i+1,:]
    dvarlimited[:,i,:] = limiter_mc(dvarup,dvardown)

    h = invdlatdist_fb[:,3:nlats-1,:]*2
    dvarup = h*(factor1*(var[:, 3:nlats-1,:]-var[:, 2:nlats-2,:])\
               -factor2*(var[:, 4:nlats,:]-var[:, 1:nlats-3,:]))
    h = invdlatdist_fb[:, 2:nlats-2,:]*2
    dvardown = h*(factor1*(var[:, 2:nlats-2,:]-var[:, 1:nlats-3,:])\
                 -factor2*(var[:, 3:nlats-1,:]-var[:, 0:nlats-4,:]))
    dvarlimited[:,1:nlats-3,:] = limiter_mc(dvarup,dvardown)

    i = nlats-3
    dvarup = (var[:, i+2,:]-var[:, i+1,:])*invdlatdist_fb[:, i+2,:]
    h = invdlatdist_fb[:,i+1,:]*2
    dvardown = h*(factor1*(var[:,i+1,:]-var[:, i,:])-factor2*(var[:,i+2,:]-var[:, i-1,:]))
    dvarlimited[:,i,:] = limiter_mc(dvarup,dvardown)

    varleft = var[:,1:nlats-2,:]+0.5*dvarlimited[:, :-1,:]*dlatdist_fb[:,2:nlats-1,:]
    varright= var[:,2:nlats-1,:]-0.5*dvarlimited[:, 1:,:]*dlatdist_fb[:,2:nlats-1,:]

    gradvar[:,2:-2,:] = \
            0.5\
            * ( varleft[:, 1:,:]  + varright[:, 1:,:]\
            - varleft[:, :-1,:] - varright[:, :-1,:] )\
            * invdlatdist_gb[:, 2:-2,:]

    return gradvar


def calc_rusanov_lons(lons, lats, alts, var):
    """
    lons: two dimensional with ghost cells (lon, lat)
    alts: two dimensional with ghost cells (lon, lat)
    var: two dimensional with ghost cells (lon, lat)
    """
    nlons, nlats, nalts = (var.shape[k] for k in [0,1,2])
    dlondist_gb = np.ones(var.shape)
    dlondist_fb = np.ones(var.shape)
    dvarlimited = np.ones([nlons-2,nlats,nalts])
    varleft = np.ones([nlons-3,nlats,nalts])
    varright = np.ones([nlons-3,nlats,nalts])
    gradvar = np.ones(var.shape)

    Re = 6371*1000
    rdist = Re+alts
    difflon = np.diff(lons, axis=0)
    coslats = np.cos(lats)
    coslats[coslats<0.01]=0.01


    dlondist_gb[1:-1, 1:-1,:] = \
            0.5*(difflon[1:,1:-1,:] + difflon[:-1,1:-1,:])\
            * (rdist[1:-1, 1:-1,:]) * coslats[1:-1, 1:-1,:]
    dlondist_gb[:, 0,:] = dlondist_gb[:,1,:]
    dlondist_gb[:,-1,:] = dlondist_gb[:,-2,:]
    dlondist_gb[0, :,:] = dlondist_gb[1,:,:]
    dlondist_gb[-1,:,:] = dlondist_gb[-2,:,:]

    dlondist_fb[1:-1, 1:-1,:] = \
            (difflon[:-1, 1:-1,:]) \
            * 0.5 * (rdist[1:-1, 1:-1,:]+ rdist[0:-2, 1:-1,:])\
            * coslats[1:-1, 1:-1,:]
    dlondist_fb[:, 0,:] = dlondist_fb[:,1,:]
    dlondist_fb[:,-1,:] = dlondist_fb[:,-2,:]
    dlondist_fb[0, :,:] = dlondist_fb[1,:,:]
    dlondist_fb[-1,:,:] = dlondist_fb[-2,:,:]

    invdlondist_gb = 1.0/dlondist_gb
    invdlondist_fb = 1.0/dlondist_fb

    factor1 = 0.6250000
    factor2 = 0.0416667

    i = 0
    h = invdlondist_fb[i+2,:,:]*2
    dvarup = h*(factor1*(var[i+2,:,:]-var[i+1,:,:]) - factor2*(var[i+3,:,:]-var[i,:,:]))
    dvardown = (var[i+1,:,:] -var[i,:,:])*invdlondist_fb[i+1,:,:]
    dvarlimited[i,:,:] = limiter_mc(dvarup,dvardown)

    h = invdlondist_fb[3:nlons-1,:,:]*2
    dvarup = h*(factor1*(var[3:nlons-1,:,:]-var[2:nlons-2,:,:])\
               -factor2*(var[4:nlons,:,:]-var[1:nlons-3,:,:]))
    h = invdlondist_fb[2:nlons-2,:,:]*2
    dvardown = h*(factor1*(var[2:nlons-2,:,:]-var[1:nlons-3,:,:])\
                 -factor2*(var[3:nlons-1,:,:]-var[0:nlons-4,:,:]))
    dvarlimited[1:nlons-3,:,:] = limiter_mc(dvarup,dvardown)

    i = nlons-3
    dvarup = (var[i+2,:,:]-var[i+1,:,:])*invdlondist_fb[i+2,:,:]
    h = invdlondist_fb[i+1,:,:]*2
    dvardown = h*(factor1*(var[i+1,:,:]-var[i,:,:])-factor2*(var[i+2,:,:]-var[i-1,:,:]))
    dvarlimited[i,:,:] = limiter_mc(dvarup,dvardown)

    varleft = var[1:nlons-2,:,:]+0.5*dvarlimited[:-1,:,:]*dlondist_fb[2:nlons-1,:,:]
    varright= var[2:nlons-1,:,:]-0.5*dvarlimited[1:,:,:]*dlondist_fb[2:nlons-1,:,:]

    gradvar[2:-2,:,:] = \
            0.5\
            * ( varleft[1:,:,:]  + varright[1:,:,:]\
            - varleft[:-1,:,:] - varright[:-1,:,:] )\
            * invdlondist_gb[2:-2, :,:]
    return gradvar


def calc_rusanov_alts_ausm(alts, var):
    factor1, factor2 = 0.625, 0.0416667
    nlons, nlats, nalts = (var.shape[k] for k in [0,1,2])

    dvarlimited = np.ones([nlons,nlats,nalts-2])*np.nan
    varleft = np.ones([nlons,nlats,nalts-3])*np.nan
    varright = np.ones([nlons,nlats,nalts-3])*np.nan
    gradvar = np.ones([nlons,nlats,nalts])*np.nan

    dalt_f = np.ones(var.shape)*np.nan
    dalt_gb = np.ones(var.shape)*np.nan

    diffalts = np.diff(alts,axis=2)
    dalt_f[:,:,1:] = diffalts
    dalt_f[:,:,0] = dalt_f[:,:,1]
    invdalt_f = 1.0/dalt_f
    dalt_gb[:,:,1:-1] = 0.5*(diffalts[:,:,1:]+diffalts[:,:,:-1])
    dalt_gb[:,:,0] = dalt_gb[:,:,1]
    dalt_gb[:,:,-1] = dalt_gb[:,:,-2]

    i = 0
    dvarup = (var[:,:,i+2]-var[:,:,i+1])*invdalt_f[:,:,i+2]
    dvardown = (var[:,:,i+1]-var[:,:,i])*invdalt_f[:,:,i+1]
    dvarlimited[:,:,i] = limiter_mc(dvarup,dvardown)

    h = invdalt_f[:,:,3:nalts-1]*2
    dvarup = h*(factor1*(var[:,:,3:nalts-1] - var[:,:,2:nalts-2])\
               -factor2*(var[:,:,4:nalts] - var[:,:,1:nalts-3]))
    h = invdalt_f[:,:,2:nalts-2]*2
    dvardown = h*(factor1*(var[:,:,2:nalts-2] - var[:,:,1:nalts-3])\
                 -factor2*(var[:,:,3:nalts-1] - var[:,:,0:nalts-4]))
    dvarlimited[:,:,1:nalts-3] = limiter_mc(dvarup,dvardown)

    i = nalts-3
    dvarup = (var[:,:,i+2]-var[:,:,i+1])*invdalt_f[:,:,i+2]
    dvardown = (var[:,:,i+1]-var[:,:,i])*invdalt_f[:,:,i+1]
    dvarlimited[:,:,i] = limiter_mc(dvarup,dvardown)

    varleft = var[:,:,1:nalts-2]+0.5*dvarlimited[:,:,:-1]*dalt_f[:,:,2:nalts-1]
    varright = var[:,:,2:nalts-1]-0.5*dvarlimited[:,:,1:]*dalt_f[:,:,2:nalts-1]

    gradvar[:,:,2:-2] = 0.5 * \
            (varleft[:,:,1:]+varright[:,:,1:]\
            -varleft[:,:,:-1]-varright[:,:,:-1])\
            /dalt_gb[:,:,2:-2]
    return gradvar


def calc_div_hozt(lons, lats, alts, varn, vare):
    Re = 6371*1000
    rdist = Re+alts
    tanlats = np.tan(lats)
    tanlats[tanlats>100]=100
    divvar = calc_rusanov_lats(lats, alts, varn)\
           + calc_rusanov_lons(lons,lats,alts,vare)\
           - tanlats*varn/rdist
    return divvar


def calc_div_vert(alts, var):
    Re = 6371*1000
    rdist = Re+alts
    divvar = calc_rusanov_alts_ausm(alts, var)\
           + 2*var/rdist
    return divvar


if __name__ == '__main__':
    import gitm
    import numpy as np
    import gitm_create_coordinate as gcc
    from cartopy.util import add_cyclic_point
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    alt=400
    g1 = gitm.GitmBin('/home/guod/simulation_output/momentum_analysis/'
            'run_shrink_iondrift_3_continue/data/3DALL_t030322_013000.bin')
    alt_ind = np.argmin(np.abs(g1['Altitude'][0,0,:]-alt*1000))
    lons, lats, alts = (g1[k] for k in ['Longitude', 'Latitude', 'Altitude'])
    veln = g1['V!Dn!N (north)']
    vele = g1['V!Dn!N (east)']+(2*np.pi)/(24*3600)*(6371*1000+alts)*np.cos(lats)
    velu = g1['V!Dn!N (up)']
    rho = g1['Rho']

    #zdata1 = velu*calc_rusanov_alts_ausm(alts,rho)
    zdata1 = calc_div_vert(alts, rho*velu)+calc_div_hozt(lons, lats, alts, rho*veln,rho*vele)
    #zdata1 = rho*calc_div_hozt(lons, lats, alts, veln, vele)
    dglat = np.array(g1['dLat'][0,2:-2,0])
    dglon = np.array(g1['dLon'][2:-2,0,0])
    zdata1, dglon = add_cyclic_point(zdata1[2:-2,2:-2,alt_ind].T, coord=dglon, axis=1)

    g2 = gitm.GitmBin('/home/guod/simulation_output/momentum_analysis/'
            'run_no_shrink_iondrift_3/data/3DALL_t030322_013000.bin')
    lons, lats, alts = (g2[k] for k in ['Longitude', 'Latitude', 'Altitude'])
    veln = g2['V!Dn!N (north)']
    vele = g2['V!Dn!N (east)']+(2*np.pi)/(24*3600)*(6371*1000+alts)*np.cos(lats)
    velu = g2['V!Dn!N (up)']
    rho = g2['Rho']
    #zdata2 = velu*calc_rusanov_alts_ausm(alts,rho)
    zdata2 = calc_div_vert(alts, rho*velu)+calc_div_hozt(lons, lats, alts, rho*veln,rho*vele)
    #zdata2 = rho*calc_div_hozt(lons, lats, alts, veln, vele)
    dglat = np.array(g2['dLat'][0,2:-2,0])
    dglon = np.array(g2['dLon'][2:-2,0,0])
    zdata2, dglon = add_cyclic_point(zdata2[2:-2,2:-2,alt_ind].T, coord=dglon, axis=1)

    dglon,dglat = np.meshgrid(dglon, dglat)
    gt = g1['time']
    centrallon = (0-(gt.hour+gt.minute/60+gt.second/3600))*15
    ax, projection = gcc.create_map(
            1, 1, 1, 'polar', nlat=-40, slat=-90, centrallon=centrallon,
            useLT=True, dlat=10)
    ax.contourf(np.array(dglon), np.array(dglat), np.array(zdata1-zdata2),
                levels=np.linspace(-5,5,21)*1e-16,transform=ccrs.PlateCarree(),
                cmap='seismic',extend='both')
    plt.show()
