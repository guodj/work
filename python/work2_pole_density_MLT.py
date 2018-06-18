#-------------------------------------------------------------------------------
# By Dongjie, USTC, Sat Sep 17 09:32:19 CST 2016
#-------------------------------------------------------------------------------

#Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.dates as mdates
import pdb   # set breakpoint
import champ_grace as cg
import omni
import os
from aacgmv2 import convert
from aacgmv2 import convert_mlt
from mpl_toolkits.axes_grid1 import ImageGrid
import datetime as dt
import cartopy.crs as ccrs
import gitm_create_coordinate as gcc



def get_sblist():
    """
    Get solar wind sector polarity reversing dates.
    Return a pd.DataFrame indexed by dates with columns sbtype, lday,
    rday and season.
    """
    sb_fname = os.environ.get('DATAPATH') + 'SBlist/SBlist.txt'
    sblist = pd.read_csv(
        sb_fname, sep='\s+', parse_dates={'dates':['year','month','day']},
        index_col='dates')
    sblist.replace(['+,-','-,+'], ['away-toward','toward-away'], inplace=True)
    doy = sblist.index.dayofyear
    sblist.loc[(doy>=35)  & (doy<=125),'season'] = 'me'
    sblist.loc[(doy>=221) & (doy<=311),'season'] = 'se'
    sblist.loc[(doy>125) & (doy<221),'season'] = 'js'
    sblist.loc[(doy>311) | (doy<35),'season'] = 'ds'
    return sblist


def get_date_polarity():
    """
    Get solar wind sector polarities and their dates.
    Return a pd.DataFrame indexed by dates with columns polarity, season.
    """

    sblist = get_sblist()
    sblist.replace(['away-toward','toward-away'], ['+,-','-,+'], inplace=True)
    date_polarity = pd.DataFrame()
    alist = []
    for row in sblist.itertuples():
        # row = ('1926-01-29', '+,-', 21, 5)
        index = (row[0] + pd.TimedeltaIndex(range(-row[2], row[3]), 'd'))
        value = list(row[2]*row[1][0]+row[3]*row[1][2])
        df = pd.DataFrame(value, index=index, columns=['polarity'])
        alist.append(df)
    date_polarity = pd.concat(alist)
    date_polarity = date_polarity.groupby(level=0).first()
    date_polarity.replace(['+','-'],['away','toward'],inplace=True)
    doy = date_polarity.index.dayofyear
    date_polarity.loc[(doy>=35)  & (doy<=125),'season'] = 'me'
    date_polarity.loc[(doy>=221) & (doy<=311),'season'] = 'se'
    date_polarity.loc[(doy>125) & (doy<221),'season'] = 'js'
    date_polarity.loc[(doy>311) | (doy<35),'season'] = 'ds'
    date_polarity.sort_index(axis=0, level=0, inplace=True)
    return date_polarity


def f3():
    def percentile(n):
        def percentile_(x):
            return np.percentile(x,n)
        percentile_.__name__ = 'percentile_%s' % n
        return percentile_
    date_polarity = get_date_polarity() # date_polarity is sorted
    date_polarity = date_polarity['2002-1-1':'2010-12-31']

    # IMF Bx, By, Bz and AE
    if False:
        print('Reading IMF data from 2002 to 2010...')
        baea = omni.get_omni(
            '2002-1-1', '2011-1-1',
            variables=['Bx', 'Bym', 'Bzm', 'AE'], res='5m')
        print('Reading finished')
        bae = [pd.DataFrame(), pd.DataFrame()]
        for k00, k0 in enumerate(['away', 'toward']):
            sbt = date_polarity[(date_polarity.polarity==k0)]
            for k11, k1 in enumerate(sbt.index):
                baet = baea[k1:(k1+pd.Timedelta('1D')-pd.Timedelta('1s'))]
                if baet.empty:
                    print('No IMF and AE data on ', k1)
                    continue
                bae[k00] = bae[k00].append(baet)
    # end of IMF data preparation

    # Grace density data.
        nsbu = np.zeros([2, 2])
        rho = [
            [pd.DataFrame(), pd.DataFrame()], [pd.DataFrame(), pd.DataFrame()]]
        for k00, k0 in enumerate(['away', 'toward']):
            sbt = date_polarity[(date_polarity.polarity==k0)]
            for k2 in sbt.index:
                rhot = cg.ChampDensity(
                    k2, k2+pd.Timedelta('1D')-pd.Timedelta('1s'),
                    satellite='grace', variables=['rho400', 'lat3'])
                if rhot.empty:
                    print('No GRACE data on ',k2)
                    continue
                for k33, k3 in enumerate([-90, 90]):  # south and north poles
                    rhott = rhot[rhot.lat3==k3].copy()
                    print([k0, k3])
                    if rhott.shape[0]<25:
                        print(
                            'There is only {:d} '
                            'data points on '.format(rhott.shape[0]), k2)
                        continue
                    rhott['rrho400'] = 100*(
                        rhott['rho400']-rhott['rho400'].mean()
                        )/rhott['rho400'].mean()
                    nsbu[k00, k33] +=1
                    rho[k00][k33] = rho[k00][k33].append(rhott)
        pd.to_pickle(
            (bae, rho, nsbu),
            os.environ.get('DATAPATH') + 'tmp/w2_f4_02.dat')
    # End of data preparation

    print('Begin figure 1')
    bdate = pd.Timestamp('2002-10-09')
    edate = pd.Timestamp('2002-10-14')
    print('Date range: ', bdate, '-->', edate)
    imf = omni.get_omni(bdate, edate, variables=['Bx', 'Bym', 'Bzm'], res='1h')
    rho = cg.ChampDensity(
        bdate, edate,
        variables=['rho400', 'lat3', 'Mlat', 'MLT'], satellite='grace')
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7.3, 6.8))
    # IMF Bx, By, Bz
    plt.sca(ax[0])
    plt.plot(imf.index, imf.Bx, 'b')
    plt.plot(imf.index, imf.Bym, 'r')
    plt.plot(imf.index, imf.Bzm, 'k')
    plt.ylim(-10, 10)
    plt.yticks(np.arange(-10, 11, 5))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.grid()
    plt.ylabel('IMF (nT)')
    plt.legend([r'$B_x$', r'$B_y$', r'$B_z$'], loc=1, ncol=3)
    plt.text(0.03, 0.87, '(a)', transform=plt.gca().transAxes)

    plt.sca(ax[1])
    rho['rho400'] /= 1e-11
    plt.plot(rho.index, rho.rho400, 'gray', lw=1)
    rhot = rho[rho.lat3==-90]
    rhott = rho[
        ((rho.Mlat<=-70) & (rho.Mlat>=-80) & (rho.MLT>=11) & (rho.MLT<=13))]
    plt.plot(rhot.index, rhot.rho400, 'k')
    plt.plot(rhott.index, rhott.rho400, 'bo', ms=5)
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator(interval=2))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.xlim(bdate, edate)
    plt.ylim(0.2, 1.4)
    plt.yticks(np.arange(0.2, 1.5, 0.2))
    dates = pd.date_range(
        bdate, edate,
        freq='1D')+pd.Timedelta('15h')+pd.Timedelta('37m')
    plt.vlines(dates, ymin=0.2, ymax=1.4, color='k', linestyle='--')
    plt.grid()
    plt.xlabel('Day of 2002')
    plt.ylabel(r'$\rho$ ($10^{-11}$ kg/m$^3$)')
    plt.text(0.03, 0.87, '(b)', transform=plt.gca().transAxes)
    plt.savefig(
        '/Users/guod/Documents/Pole_Density_MLT_Change/Figures/001_01_Case.pdf')
    print('End of figure 1\n\n')

    print('Begin figure 2')
    bae, rho, nsbu= pd.read_pickle(
        os.environ.get('DATAPATH') + 'tmp/w2_f4_02.dat')
    print('Total away and toward days [[AS, AN], [TS, TN]]: \n', nsbu)
    fig = plt.figure(figsize=(6,8))
    grid1 = ImageGrid(
        fig, [0.1, 0.4, 0.8, 0.66], [3, 2], label_mode='L',axes_pad=0.2,
        cbar_mode='edge', cbar_location='right')
    grid2 = ImageGrid(
        fig, [0.1, 0.07, 0.8, 0.44], [2, 2], label_mode='L',
        axes_pad=[0.2, 0.4], cbar_mode='edge', cbar_location='right')
    ctt = [r'$B_x$ (nT)', r'$B_y$ (nT)', r'$B_z$ (nT)']
    plb = np.array(list('abcdefghij')).reshape(2,5).T
    for k00, k0 in enumerate(['Away', 'Toward']):
        baet = bae[k00]
        baet['month'] = baet.index.month-0.5
        baet['uthour'] = baet.index.hour+0.5
        baett = baet.groupby(['month', 'uthour']).agg(np.median)
        baett = baett.reset_index()
        print('For %s polarities: ' % k0)
        for k11, k1 in enumerate(['Bx', 'Bym', 'Bzm']):
            ll = np.linspace(-3.2, 3.2, 11)
            cl = np.arange(-3, 4, 1)
            if k1=='AE':
                ll = np.linspace(0, 300, 11)
                cl = np.arange(0, 301, 100)
            #plt.sca(ax[k11, k00])
            plt.sca(grid1[k00+k11*2])
            baettt = baett.pivot('month', 'uthour', k1)
            # Extend month and ut
            baettt.loc[:, baettt.columns[0]-1] = baettt.loc[:, baettt.columns[-1]]
            baettt.loc[:, baettt.columns[-2]+1] = baettt.loc[:, baettt.columns[0]]
            baettt = baettt.sort_index(axis=1)
            baettt.loc[baettt.index[0]-1, :] = baettt.loc[baettt.index[-1], :]
            baettt.loc[baettt.index[-2]+1, :] = baettt.loc[baettt.index[0], :]
            baettt = baettt.sort_index(axis=0)
            x = baettt.columns # uthour
            y = baettt.index # month
            hc = plt.contourf(
                x, y, baettt, levels=ll, extend='neither', cmap='seismic')
            print('  Average {:s} is: {:5.1f}'.format(k1, baettt.mean().mean()))
            if k1=='Bzm':
                print('  Bz max: {:5.1f}'.format(baettt.max().max()))
                print('  Bz min: {:5.1f}'.format(baettt.min().min()))
            plt.xlim(0,24)
            plt.xticks(np.arange(0,25,6),[])
            plt.yticks(np.arange(0.5,12.5,1),['','',3,'','',6,'','',9,'','',12])
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
            plt.tick_params(axis='x', which='major', direction='out', length=4)
            plt.tick_params(axis='y', which='major', direction='out', length=0, pad=6)
            plt.tick_params(axis='x', which='minor', direction='out', length=2)
            plt.tick_params(axis='y', which='minor', direction='out', length=3, width=1.2)
            plt.text(
                0.1, 0.82, '('+plb[k11, k00]+')',bbox=dict(facecolor='grey', alpha=0.5),
                transform=plt.gca().transAxes)
            if k11 == 0:
                plt.title(k0)
            if k00 == 0:
                plt.ylabel('Month')
            plt.ylim(-0.000001,12.0000001)
            grid1.cbar_axes[k11].colorbar(hc, ticks=cl)
            grid1.cbar_axes[k11].set_ylabel(ctt[k11])
        for k11, k1 in enumerate(['south', 'north']):
            rhot = rho[k00][k11]
            rhot['month'] = rhot.index.month
            rhot['uthour'] = rhot.index.hour+0.5
            rhott = rhot.groupby(['month', 'uthour']).agg(np.median)
            rhott = rhott.reset_index()
            plt.sca(grid2[k00+k11*2])
            rhottt = rhott.pivot('month', 'uthour', 'rrho400')
            # extend month and ut
            rhottt.loc[:, rhottt.columns[0]-1] = rhottt.loc[:, rhottt.columns[-1]]
            rhottt.loc[:, rhottt.columns[-2]+1] = rhottt.loc[:, rhottt.columns[0]]
            rhottt = rhottt.sort_index(axis=1)
            rhottt.loc[rhottt.index[0]-1, :] = rhottt.loc[rhottt.index[-1], :]
            rhottt.loc[rhottt.index[-2]+1, :] = rhottt.loc[rhottt.index[0], :]
            rhottt = rhottt.sort_index(axis=0)
            hc = plt.contourf(
                x, y, rhottt, levels=np.linspace(-22, 22, 11), cmap='seismic')
            print('  ', k1, ' density max (%): ', rhottt.max().max())
            print('  ', k1, ' density min (%): ', rhottt.min().min())
            if k1 is 'south':
                #plt.axvline(15+37/60, 0, 1, c='k', ls='--')
                utx = np.arange(0,25,6)
                uts = [convert_mlt(19,dtime=dt.datetime(2003,2,23,k)) for k in utx%24]
                [plt.text(k1, 13, '%.0f'%k2, horizontalalignment='center')\
                    for k1, k2 in zip(utx, uts)]
                if k00==0:
                    plt.text(-0.2,1.08,'MLT',transform=plt.gca().transAxes)
            if k1 is 'north':
                #plt.axvline(5+25/60, 0, 1, c='k', ls='--')
                utx = np.arange(0,25,6)
                utn = [convert_mlt(170,dtime=dt.datetime(2003,2,23,k)) for k in utx%24]
                [plt.text(k1, 13, '%.0f'%k2, horizontalalignment='center')\
                    for k1, k2 in zip(utx, utn)]
                if k00==0:
                    plt.text(-0.2,1.08,'MLT',transform=plt.gca().transAxes)
            plt.xlim(0, 24)
            plt.xticks(np.arange(0, 25, 6))
            plt.ylim(-0.001, 12.001)
            plt.yticks(
                np.arange(0.5, 12.5),
                ['', '', 3, '', '', 6, '', '', 9, '', '', 12])
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
            plt.tick_params(axis='x', which='major', direction='out', length=4)
            plt.tick_params(
                axis='y', which='major', direction='out', length=0, pad=6)
            plt.tick_params(axis='x', which='minor', direction='out', length=2)
            plt.tick_params(
                axis='y', which='minor', direction='out', length=3, width=1.2)
            plt.text(
                0.1, 0.82, '('+plb[k11+3, k00]+')',
                bbox=dict(facecolor='grey', alpha=0.5),
                transform=plt.gca().transAxes)
            if k00 == 0:
                plt.ylabel('Month')
            if k11 == 1:
                plt.xlabel('UT (hour)')
            grid2.cbar_axes[k11].colorbar(hc, ticks=np.arange(-20,21,10))
            grid2.cbar_axes[k11].set_ylabel(
                (r'S' if k11==0 else 'N') +r', $\delta\rho$ (%)')
    plt.savefig(
        '/Users/guod/Documents/Pole_Density_MLT_Change/'
        'Figures/002_Statistical_Results.pdf')
    print('End of figure 2\n\n')

    print('Begin figure 3')
    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 4))
    plt.subplots_adjust(
        left=0.10, right=0.95, top=0.90, bottom=0.15, wspace=0.23, hspace=0.16)
    plb = ['(a)', '(b)']
    for k00, k0 in enumerate(['Solar Maximum', 'Solar Minimum']):
        print(k0, ':')
        rhot = rho[0][0] #  away sectors, south pole
        if 'max' in k0.lower():
            rhott = rhot['2002-1-1':'2003-12-31'].copy()
            print('  {:d} days'.format(len(np.unique(rhott.index.date))))
            tit = 'Year: 2002-2003'
        else:
            rhott = rhot['2009-1-1':'2010-12-31'].copy()
            print('  {:d} days'.format(len(np.unique(rhott.index.date))))
            tit = 'Year: 2009-2010'
        rhott = rhott[(rhott.index.month>=9) & (rhott.index.month<=10)].copy()
        rhott['uthour'] = rhott.index.hour+0.5
        rhottt = rhott.groupby(['uthour'])['rrho400'].agg(
            [np.median, percentile(25), percentile(75)])
        rhottt.columns = ['median', 'p25', 'p75']
        plt.sca(ax[k00])
        # Extend ut
        rhottt.loc[rhottt.index[0]-1, :] = rhottt.loc[rhottt.index[-1], :]
        rhottt.loc[rhottt.index[-2]+1, :] = rhottt.loc[rhottt.index[0], :]
        rhottt = rhottt.sort_index(axis=0)
        hp = plt.plot(rhottt.index, rhottt['median'], 'b')
        print('  Density max (%): ', rhottt['median'].max())
        print('  Density min (%): ', rhottt['median'].min())
        plt.plot(
            rhottt.index, rhottt.p25,'gray', rhottt.index, rhottt.p75,
            'gray', linestyle='--')
        plt.grid()
        plt.xlim(0, 24)
        plt.xticks(np.arange(0, 25, 6), [])
        utx = np.arange(0,25,6)
        uts = [convert_mlt(19,dtime=dt.datetime(2003,2,23,k)) for k in utx%24]
        plt.text(-3.5, -34, 'UT')
        [plt.text(k1, -34, '%.0f'%k1, horizontalalignment='center')
            for k1 in utx]
        plt.text(-3.5, -38, 'MLT')
        [plt.text(k1, -38, '%.0f'%k2, horizontalalignment='center')
            for k1, k2 in zip(utx, uts)]

        plt.ylim(-30, 30)
        plt.yticks(np.arange(-30, 31, 10))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.tick_params(axis='both', which='major')
        plt.tick_params(axis='both', which='minor')
        if k00 == 0:
            plt.ylabel(r'South, $\delta\rho$ (%)')
        plt.title(tit)
        plt.text(0.03, 0.91, plb[k00], transform=plt.gca().transAxes)
    print('End of figure 3\n\n')
    plt.savefig(
        '/Users/guod/Documents/Pole_Density_MLT_Change/'
        'Figures/003_Solar_Activity_Dependence.pdf')

    print('Begin figure 4')
    bdate = pd.Timestamp('2002-10-09')
    edate = pd.Timestamp('2002-10-14')
    print('Date range: ', bdate, '-->', edate)
    rho = cg.ChampDensity(
        bdate, edate,
        variables=['rho400', 'lat3', 'Mlat', 'MLT'], satellite='grace')
    rho['Dist_Cusp'] = 6371*np.arccos(
        np.sin(rho['Mlat']/180*np.pi)*np.sin(-75/180*np.pi) +
        np.cos(rho['Mlat']/180*np.pi)*np.cos(-75/180*np.pi) *
            np.cos((rho['MLT']-12)/12*np.pi))
    fig, ax = plt.subplots(1, 1, figsize=(7.3, 4.8))
    rho['rho400'] /= 1e-11
    #plt.plot(rho.index, rho.rho400, 'gray', lw=1)
    rhot = rho[rho.lat3==-90]
    plt.plot(rhot['Dist_Cusp'], rhot['rho400'], '.k')
    plt.xlim(0, 3500)
    plt.ylim(0.2, 1.4)
    plt.yticks(np.arange(0.2, 1.5, 0.2))
    plt.xticks(np.arange(0,4000,500))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.grid()
    plt.xlabel('Pole-Cusp Distance (km)')
    plt.ylabel(r'$\rho$ ($10^{-11}$ kg/m$^3$)')
    print('End of figure 4\n\n')
    plt.savefig(
        '/Users/guod/Documents/Pole_Density_MLT_Change/'
        'Figures/001_02_Dist_Cusp.pdf')
    return

def f4():
    plt.figure(figsize=(6.88,6.74))
    #geographic coordinates
    ax1, projection1 = gcc.create_map(
        2, 2, 1, 'pol', 90, 50, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
    ax1.plot([45,45],[50,90], 'k--',transform=ccrs.PlateCarree())
    ax1.plot([225,225],[50,90], 'k--',transform=ccrs.PlateCarree())
    ax1.plot([105,105],[50,90], 'k--',transform=ccrs.PlateCarree())
    ax1.plot([285,285],[50,90], 'k--',transform=ccrs.PlateCarree())
    ax1.scatter(
        0, 90, color='r', transform=ccrs.PlateCarree(), zorder=10,
        label='North Pole')
    ax1.text(0,1,'(a)', transform = plt.gca().transAxes)
    plt.legend(loc=[0.5,1.1])

    ax2, projection2 = gcc.create_map(
        2, 2, 2, 'pol', -50, -90, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
    ax2.scatter(
        0, -90, color='b', transform=ccrs.PlateCarree(),label='South Pole')
    ax2.text(0,1,'(b)', transform = plt.gca().transAxes)
    plt.legend(loc=[-0.1,1.1])

    #geomagnetic coordinates
    ax3, projection3 = gcc.create_map(
        2, 2, 3, 'pol', 90, 50, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
    mlatn,mlonn = convert(90,0,0,date=dt.date(2002,3,21))
    for k in range(24):
        mltn = convert_mlt(mlonn[0],dtime=dt.datetime(2003,3,21,k))
        ax3.scatter(mltn*15,mlatn[0],color='r',transform=ccrs.PlateCarree())
    ax3.scatter(180,75,s=50,c='k',marker='x',transform=ccrs.PlateCarree())
    ax3.text(0,1,'(c)', transform = plt.gca().transAxes)

    ax4, projection4 = gcc.create_map(
        2, 2, 4, 'pol', -50, -90, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
    mlats,mlons = convert(-90,0,0,date=dt.date(2002,3,21))
    for k in range(24):
        mlts = convert_mlt(mlons[0],dtime=dt.datetime(2003,3,21,k))
        ax4.scatter(mlts*15,mlats[0],color='b',transform=ccrs.PlateCarree())
    ax4.scatter(180,-75,s=50,c='k',marker='x',transform=ccrs.PlateCarree())
    ax4.text(0,1,'(d)', transform = plt.gca().transAxes)
    plt.savefig(
        '/Users/guod/Documents/Pole_Density_MLT_Change/Figures/'
        '000_Pole_Feature.pdf')
    return
# END
#-------------------------------------------------------------------------------
if __name__=='__main__':
    plt.close('all')
    a = f3()
    #f4()
    plt.show()
