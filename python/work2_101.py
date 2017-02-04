#-------------------------------------------------------------------------------
# For the second project
# By Dongjie, USTC, Sat Sep 17 09:32:19 CST 2016
#-------------------------------------------------------------------------------

#Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pdb   # set breakpoint
import champ_grace as cg
import omni
import os


def get_sblist():
    """ Get solar wind sector polarity reversing dates

    Args:
        No input
    Returns:
        A dataframe indexed by dates with columns sbtype, lday, rday and season.
    """
    sb_fname = os.environ.get('DATAPATH') + 'SBlist/SBlist.txt'
    sblist = pd.read_csv(
        sb_fname,
        sep='\s+',
        parse_dates={'dates':['year','month','day']},
        index_col='dates')
    sblist.replace(['+,-','-,+'], ['away-toward','toward-away'], inplace=True)
    doy = sblist.index.dayofyear
    sblist.ix[(doy>=35)  & (doy<=125),'season'] = 'me'
    sblist.ix[(doy>=221) & (doy<=311),'season'] = 'se'
    sblist.ix[(doy>125) & (doy<221),'season'] = 'js'
    sblist.ix[(doy>311) | (doy<35),'season'] = 'ds'
    return sblist


def get_date_polarity():
    """ Get solar wind sector polarities and their dates.

    Args:
        no input
    Returns:
        date_polarity: pd.DataFrame with columns 'polarity', 'season'.
            The index is dates
    """

    sblist = get_sblist()
    sblist.replace(['away-toward','toward-away'], ['+,-','-,+'], inplace=True)
    date_polarity = pd.DataFrame()
    alist = []
    for row in sblist.itertuples():
        # row = ('1926-01-29', 'away-toward', 21, 5)
        index = (row[0] + pd.TimedeltaIndex(range(-row[2], row[3]), 'd'))
        value = list(row[2]*row[1][0]+row[3]*row[1][2])
        df = pd.DataFrame(value, index=index, columns=['polarity'])
        alist.append(df)
    date_polarity = pd.concat(alist)
    date_polarity = date_polarity.groupby(level=0).first()
    date_polarity.replace(['+','-'],['away','toward'],inplace=True)
    doy = date_polarity.index.dayofyear
    date_polarity.ix[(doy>=35)  & (doy<=125),'season'] = 'me'
    date_polarity.ix[(doy>=221) & (doy<=311),'season'] = 'se'
    date_polarity.ix[(doy>125) & (doy<221),'season'] = 'js'
    date_polarity.ix[(doy>311) | (doy<35),'season'] = 'ds'
    date_polarity.sort_index(axis=0, level=0, inplace=True)
    return date_polarity


def f1():
    '''Imf, AE and density variations during epoch days of sblist
    Before running f1, run f2 firstly.
    '''
    sblist = get_sblist()
    sblist = sblist['2002-8-1':'2010-6-27']
    index_name = ['Bx', 'Bym', 'Bzm', 'AE']
    tl = ['$B_x$ (nT)','$B_y$ (nT)','$B_z$ (nT)','AE']
    levels = [np.linspace(-4, 4, 11), np.linspace(-4, 4, 11),
              np.linspace(-4, 4, 11), np.linspace(0, 400, 11)]
    cticks = [np.arange(-4, 5, 2), np.arange(-4, 5, 2),
              np.arange(-4, 5, 2), np.arange(0, 401, 100)]
    data = [pd.DataFrame(),pd.DataFrame()]
    nn = [0, 0]
    if False:  # data preparation
        for k00, k0 in enumerate(['away-toward', 'toward-away']):
            sbtmp = sblist[sblist.sbtype==k0]
            for k1 in sbtmp.index:
                bdate = k1-pd.Timedelta('3D')
                edate = k1+pd.Timedelta('3D')
                data_tmp = get_omni(bdate, edate, index_name, res='1h')
                if not data_tmp.empty:
                    nn[k00] = nn[k00]+1
                    print(nn, k1, k0)
                    data_tmp['epochday'] = (
                            data_tmp.index - k1)/pd.Timedelta('1D')
                    data[k00] = data[k00].append(data_tmp)
        pd.to_pickle(data, os.environ.get('DATAPATH') + 'tmp/w2_07.dat')
    # END of data preperation
    fig,ax = plt.subplots(6,2,sharey=True,figsize=(7.3,9))
    # IMF and AE index
    data = pd.read_pickle(os.environ.get('DATAPATH') + 'tmp/w2_07.dat')
    datagroup = [
            data[k].groupby([data[k].index.month,np.floor(data[k].epochday*24)])
            for k in [0,1]]
    datagroup = [datagroup[k].median() for k in [0,1]]
    for k in [0,1]:
        datagroup[k].index.names = ('month', 'epochhour')
        datagroup[k] = datagroup[k].reset_index().pivot(
                index='epochhour', columns='month')
        for k1 in datagroup[k].columns.levels[0]:
            datagroup[k][k1, 13] = datagroup[k][k1, 1]
    for k00, k0 in enumerate(['Away-Toward', 'Toward-Away']):
        for k11,k1 in enumerate(index_name):
            plt.sca(ax[k11,k00])
            data1 = datagroup[k00][k1]
            hc2 = plt.contourf(data1.columns, data1.index/24, data1.values,
                    levels=levels[k11], cmap='bwr', extend='both')
            plt.xlim([1,13])
            plt.xticks(np.arange(1, 14))
            plt.gca().set_xticklabels('')
            plt.gca().set_xticks(np.arange(1.5, 13.5, 1),  minor=True)
            plt.ylim([-3, 3])
            plt.yticks(np.arange(-3, 4, 1), fontsize=13)
            plt.tick_params(axis='both', which='major',
                            direction='out', length=3.5)
            plt.tick_params(axis='both', which='minor', length=0)
            if k00 is 1:
                axpo = np.array(plt.gca().get_position())
                cax = plt.gcf().add_axes(
                        (axpo[1,0]-0.005,axpo[0,1],0.01,axpo[1,1]-axpo[0,1]))
                cbar = plt.colorbar(mappable=hc2,cax=cax,ticks=cticks[k11])
                cbar.set_label(tl[k11])
                plt.tick_params('both', length=0)
    # Density
    data = pd.read_pickle(os.environ.get('DATAPATH') + 'tmp/w2_13.dat')
    datagroup = [[None, None], [None, None]]
    rhoname = 'rrho400'
    tl = [r'N $\rho_r$ (%)', r'S $\rho_r$ (%)']
    levels = [np.linspace(-45, 45, 11), np.linspace(-45, 45, 11)]
    cticks = [np.arange(-45, 46, 15), np.arange(-45, 46, 15)]
    for k0 in range(2):
        for k1 in range(2):
            datagrouptmp = data[k0][k1].groupby(
                    [data[k0][k1].index.month,
                     np.floor(data[k0][k1].epochday*48)])
            datagroup[k0][k1] = datagrouptmp.median()
            datagroup[k0][k1].index.names = ('month', 'epochhalfhour')
            datagroup[k0][k1] = datagroup[k0][k1].reset_index().pivot(
                    index='epochhalfhour', columns='month')
            for k2 in datagroup[k0][k1].columns.levels[0]:
                datagroup[k0][k1][k2, 13] = datagroup[k0][k1][k2, 1]
    for k00, k0 in enumerate(['Away-Toward', 'Toward-Away']):
        for k11, k1 in enumerate(['N', 'S']):
            plt.sca(ax[k11+len(index_name),k00])
            data1 = datagroup[k00][k11][rhoname]
            hc2 = plt.contourf(data1.columns, data1.index/48, data1.values,
                    levels=levels[k11], cmap='bwr', extend='both')
            plt.xlim([1,13])
            plt.xticks(np.arange(1, 14))
            plt.gca().set_xticklabels('')
            plt.gca().set_xticks(np.arange(1.5, 13.5, 1),  minor=True)
            plt.ylim([-3, 3])
            plt.yticks(np.arange(-3, 4, 1), fontsize=13)
            plt.tick_params(axis='both', which='major',
                            direction='out', length=3.5)
            plt.tick_params(axis='both', which='minor', length=0)
            if k00 is 1:
                axpo = np.array(plt.gca().get_position())
                cax = plt.gcf().add_axes(
                        (axpo[1,0]-0.005,axpo[0,1],0.01,axpo[1,1]-axpo[0,1]))
                cbar = plt.colorbar(mappable=hc2,cax=cax,ticks=cticks[k11])
                cbar.set_label(tl[k11])
                plt.tick_params('both', length=0)
    # Set ax
    ax[-1, 0].set_xticklabels(np.arange(1, 13), minor=True, fontsize=13)
    ax[-1, 1].set_xticklabels(np.arange(1, 13), minor=True, fontsize=13)
    title1 = ax[0,0].set_title('Away-Toward')
    title2 = ax[0,1].set_title('Toward-Away')
    title1.set_position((0.5,1.05))
    title2.set_position((0.5,1.05))
    ax[-1,0].set_xlabel('Month',fontsize=14)
    ax[-1,1].set_xlabel('Month',fontsize=14)
    plt.text(0.03,0.5,'Epoch Time (day)',fontsize=14,
             verticalalignment='center',
             transform=plt.gcf().transFigure,
             rotation='vertical')
    plt.subplots_adjust(right=0.87,wspace=0.08)
    return


def f2():
    """Near the geodetic poles, the longitude and LT are not important;
    so it is a good for research of UT variation
    """
    def percentile(n):
        def percentile_(x):
            return np.percentile(x,n)
        percentile_.__name__ = 'percentile_%s' % n
        return percentile_

    sblist = get_sblist()
    sblist = sblist['2001-1-1':'2010-12-31']
    # [[ATN,ATS],[TAN,TAS]]
    density = [[pd.DataFrame(),pd.DataFrame()],[pd.DataFrame(),pd.DataFrame()]]
    sbn = [0,0]
    if False:
        for k00,k in enumerate(['away-toward','toward-away']):
            sbtmp = sblist[sblist.sbtype==k]
            for k1 in sbtmp.index:
                #for k2 in ['champ','grace']:
                for k2 in ['grace']:  # only consider the grace
                    rho = cg.ChampDensity(
                            k1-pd.Timedelta('3D'), k1+pd.Timedelta('3D'), k2)
                    if rho.empty:
                        print('no data around',k1)
                        continue
                    # The geomagnetic activity is not considered, because cusp
                    # density enhancement occurs in both geomagnetic quiet and
                    # active conditions
                    """
                    rho = rho.add_index()
                    l1 = len(rho)
                    rho = rho[rho.Kp<=40]
                    l2 = len(rho)
                    if l2<0.8*l1:
                        print('active geomagnetic condition around', k1)
                        continue
                    """
                    #rho1 = rho[rho.lat3>=87]  # north pole
                    rho1 = rho[rho.lat3==90].copy() # only consider the grace
                    # Make sure that there is data in each day
                    if len(np.unique(rho1.index.dayofyear))!=6:
                        print('there is data gap around', k1)
                        continue
                    rho1['epochday'] = (rho1.index-k1)/pd.Timedelta('1D')
                    rho1['rrho400'] = 100*(
                            rho1.rho400-rho1['rho400'].mean()
                            )/rho1['rho400'].mean()
                    rho1['rrho'] = 100*(
                            rho1.rho-rho1['rho'].mean())/rho1['rho'].mean()
                    rho1['rmsis'] = 100*(
                            rho1.msis_rho-rho1['msis_rho'].mean()
                            )/rho1['msis_rho'].mean()
                    density[k00][0] = density[k00][0].append(
                            rho1[['epochday', 'MLT', 'Mlat',
                                  'rho', 'rrho', 'rho400', 'rrho400',
                                  'msis_rho', 'rmsis']])

                    #rho2 = rho[rho.lat3<=-87]  # south pole
                    rho2 = rho[rho.lat3==-90].copy()  # only consider the grace
                    if len(np.unique(rho2.index.dayofyear))!=6:
                        print('there is data gap around', k1)
                        continue
                    rho2['epochday'] = (rho2.index-k1)/pd.Timedelta('1D')
                    rho2['rrho400'] = 100*(
                            rho2.rho400-rho2['rho400'].mean()
                            )/rho2['rho400'].mean()
                    rho2['rrho'] = 100*(
                            rho2.rho-rho2['rho'].mean())/rho2['rho'].mean()
                    rho2['rmsis'] = 100*(
                            rho2.msis_rho-rho2['msis_rho'].mean()
                            )/rho2['msis_rho'].mean()
                    density[k00][1] = density[k00][1].append(
                            rho2[['epochday', 'MLT', 'Mlat',
                                  'rho', 'rrho', 'rho400', 'rrho400',
                                  'msis_rho', 'rmsis']])
                    sbn[k00] = sbn[k00]+1
                    print(sbn)
        pd.to_pickle(density, os.environ.get('DATAPATH') + 'tmp/w2_13.dat')
    # END of data preperation

    # Pole density variation as a function of epoch time at different seasons
    # and sbtype.
    density = pd.read_pickle(os.environ.get('DATAPATH') + 'tmp/w2_13.dat')
    fig,ax = plt.subplots(4,4,sharex=True,sharey=True,figsize=(8,8))
    # fl = [['(a1)','(a2)','(a3)','(a4)'],['(b1)','(b2)','(b3)','(b4)'],
    #       ['(c1)','(c2)','(c3)','(c4)'],['(d1)','(d2)','(d3)','(d4)']]
    fl = ['(a)', '(b)', '(c)', '(d)']
    # case number in each season catagary
    nn = np.zeros([4,4])*np.nan
    for k00,k in enumerate(['away-toward','toward-away']):
        for k11, k1 in enumerate(['N','S']):
            density1 = density[k00][k11]
            # density1 = density1['2002-1-1':'2004-12-31']
            if density1.empty:
                continue
            for k22, k2 in enumerate(['me','se','js','ds']):
                plt.sca(ax[k22,k00*2+k11])
                if k2 is 'me':
                    fp = (density1.index.month>=2) & (density1.index.month<=4)
                if k2 is 'se':
                    fp = (density1.index.month>=8) & (density1.index.month<=10)
                if k2 is 'js':
                    fp = (density1.index.month>=5) & (density1.index.month<=7)
                if k2 is 'ds':
                    fp = (density1.index.month>=11) | (density1.index.month<=1)
                density2 = density1[fp].copy()
                nn[k22, k00*2+k11] = len(np.unique(
                        density2[(density2.epochday>=0) &
                        (density2.epochday<1)].index.date))
                density2['epochbin'] = density2.epochday*24//1.5*1.5+0.75
                density2 = density2.groupby('epochbin')['rrho400'].agg(
                        [np.median, percentile(25),percentile(75)])
                density2.columns = ['median', 'p25', 'p75']
                plt.plot(density2.index/24, density2['p25'],'gray',
                         density2.index/24, density2['p75'],'gray',
                         linestyle='--',dashes=(2,1),linewidth=1)
                plt.plot(density2.index/24, density2['median'],'b',linewidth=2)
                plt.xlim(-3,3)
                plt.xticks(np.arange(-3,4,1))
                #plt.gca().xaxis.set_minor_locator(AutoMinorLocator(4))
                if k1 is 'S':
                    plt.vlines(np.arange(-3, 3)+15.5/24, -30, 60,
                            linestyle='--', linewidth=1, color='r')
                if k1 is 'N':
                    plt.vlines(np.arange(-3, 3)+5.5/24, -30, 60,
                            linestyle='--', linewidth=1,color='r')
                plt.ylim(-30,60)
                plt.yticks(np.arange(-30,61,30))
                #plt.grid(which='minor',dashes=(4,1))
                #plt.grid(which='major',axis='y',dashes=(4,1))
                plt.grid(axis='y',dashes=(4,1))
                plt.tick_params(
                        axis='both',which='major',direction='out',length=4)
                if k00*2+k11==0:
                    plt.ylabel(r'$\rho_r$ (%)')
                if k22==3:
                    plt.xlabel('Epoch Time (day)',fontsize=12)
                plt.text(0.1,0.8,k1,transform=plt.gca().transAxes)
                if k22==0:
                    plt.text(
                            0,1.07,fl[k00*2+k11], transform=plt.gca().transAxes)
    plt.subplots_adjust(left=0.11,wspace=0.11, hspace=0.12)
    plt.text(0.21,0.95,'Away - Toward',transform=plt.gcf().transFigure)
    plt.text(0.61,0.95,'Toward - Away',transform=plt.gcf().transFigure)
    plt.text(0.91,0.8,'Feb - Apr',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.59,'Aug - Oct',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.38,'May - Jul',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.17,'Nov - Jan',transform=plt.gcf().transFigure,fontsize=11)
    print(nn)

    # Density, By and Bz variations versus UT
    fig, ax = plt.subplots(3, 2, sharex=True, sharey='row', figsize=(6.35, 7.02))
    density1 = density[1][1]  # Away-Toward, south pole
    # December solstice
    density1 = density1[(density1.index.month>=11) | (density1.index.month<=1)]
    fl = ['(a)', '(b)']
    for k00,k0 in enumerate(['Away', 'Toward']):
        plt.sca(ax[0, k00])
        if k0 is 'Away':
            density2 = density1[density1.epochday<=0].copy()
        if k0 is 'Toward':
            density2 = density1[density1.epochday>0].copy()
        density2['epochhourbin'] = (density2.epochday % 1)*24//1.5*1.5+0.75
        density2 = density2.groupby('epochhourbin')['rrho400'].agg(
                [np.median, percentile(25),percentile(75)])
        density2.columns = ['median', 'p25', 'p75']
        density2.loc[-0.75] = density2.loc[23.25]
        density2.loc[24.75] = density2.loc[0.75]
        density2.sort_index(inplace=True)
        plt.plot(density2.index, density2['p25'],'gray',
                 density2.index, density2['p75'],'gray',
                 linestyle='--',dashes=(2,1),linewidth=1)
        plt.plot(density2.index, density2['median'], 'b', linewidth=2)
        plt.xlim(0, 24)
        plt.xticks(np.arange(0, 25, 4))
        plt.axvline(15.5, linestyle='--', linewidth=1, color='r')
        plt.ylim(-30,30)
        plt.yticks(np.arange(-30,31,15))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(4))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(3))
        plt.tick_params(axis='both', which='major', length=4)
        plt.grid(dashes=(4,1))
        plt.title(k0)
        plt.text(0,1.05,fl[k00], transform=plt.gca().transAxes)
    # IMF By, Bz and AE
    data = pd.read_pickle(os.environ.get('DATAPATH') + 'tmp/w2_07.dat')
    data = data[0]  # away-toward
    # December solstice
    data = data[(data.index.month>=11) | (data.index.month<=1)]
    fl = [['(c)', '(d)'], ['(e)', '(f)']]
    for k00,k0 in enumerate(['Away', 'Toward']):
        for k11, k1 in enumerate(['Bym', 'Bzm']):
            plt.sca(ax[k11+1, k00])
            if k0 is 'Away':
                data1 = data[data.epochday<=0].copy()
            if k0 is 'Toward':
                data1 = data[data.epochday>0].copy()
            data1['epochhourbin'] = (data1.epochday % 1)*24//1*1+0.5
            data1 = data1.groupby('epochhourbin')[k1].agg(
                    [np.median, percentile(25), percentile(75)])
            data1.columns = ['median', 'p25', 'p75']
            data1.loc[-0.5] = data1.loc[23.5]
            data1.loc[24.5] = data1.loc[0.5]
            data1.sort_index(inplace=True)
            plt.plot(data1.index, data1['p25'],'gray',
                     data1.index, data1['p75'],'gray',
                     linestyle='--',dashes=(2,1),linewidth=1)
            plt.plot(data1.index, data1['median'], 'b', linewidth=2)
            plt.xlim(0, 24)
            plt.xticks(np.arange(0, 25, 4))
            plt.ylim(-4,4)
            plt.yticks(np.arange(-4,5,2))
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(4))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
            plt.tick_params(axis='both', which='major', length=4)
            plt.grid(dashes=(4,1))
            plt.text(0,1.05,fl[k11][k00], transform=plt.gca().transAxes)
    ax[0,0].set_ylabel(r'$\rho_r$ (%)',fontsize=14)
    ax[1,0].set_ylabel(r'$B_y$ (nT)',fontsize=14)
    ax[2,0].set_ylabel(r'$B_z$ (nT)',fontsize=14)
    [ax[-1, k].set_xlabel('UT (hour)',fontsize=14) for k in range(2)]
    plt.subplots_adjust(left=0.15,wspace=0.08,hspace=0.24,bottom=0.1)

    # Density variations at solar maximum and minimum.
    fig,ax = plt.subplots(2,1,sharex=True,sharey=True,figsize=(6.35,7.02))
    density1 = density[0][1] # for away-toward and south pole
    fl = ['(a)','(b)']
    nn = [0,0]
    for k00,k0 in enumerate(['Solar maximum','Solar minimum']):
        plt.sca(ax[k00])
        if k0 is 'Solar maximum':
            density2 = density1['2002-1-1':'2004-12-31']
        if k0 is 'Solar minimum':
            density2 = density1['2005-1-1':'2010-12-31']
        density2 = density2[
                (density2.index.month>=8) & (density2.index.month<=10)]
        nn[k00] = len(np.unique(
                density2[(density2.epochday>=0) &
                         (density2.epochday<1)].index.date))
        density2['epochbin'] = density2.epochday*24//1.5*1.5+0.75
        density2 = density2.groupby('epochbin')['rrho400'].agg(
                [np.median, percentile(25),percentile(75)])
        density2.columns = ['median', 'p25', 'p75']
        plt.plot(density2.index/24, density2['p25'],'gray',
                 density2.index/24, density2['p75'],'gray',
                 linestyle='--',dashes=(2,1),linewidth=1)
        plt.plot(density2.index/24, density2['median'],'b',linewidth=2)
        plt.xlim(-3,3)
        plt.xticks(np.arange(-3,4,1))
        plt.vlines(np.arange(-3, 3)+15.5/24, -30, 60,
                   linestyle='--', linewidth=1, color='r')
        plt.ylim(-30,60)
        plt.yticks(np.arange(-30,61,30))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(3))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.tick_params(axis='both', which='major', length=5)
        #plt.grid(which='minor',dashes=(4,1))
        #plt.grid(which='major',axis='y',dashes=(4,1))
        plt.grid(axis='y',dashes=(4,1))
        plt.ylabel(r'$\Delta\rho$ (%)',fontsize=14)
        if k00==1:
            plt.xlabel('Epoch Time (day)',fontsize=14)
        if k00==0:
            a = plt.title('Year: 02 - 04')
        if k00==1:
            a = plt.title('Year: 05 - 10')
        a.set_position((0.5,1.06))
        plt.text(0.1,0.8,'S',transform=plt.gca().transAxes)
        plt.text(0,1.05,fl[k00], transform=plt.gca().transAxes)
    print(nn)
    plt.subplots_adjust(left=0.15,wspace=0.04,hspace=0.24,bottom=0.1)

    # Magnetic local time changes at two poles as a function of UT
    fig,ax = plt.subplots(2,1,sharex=True,sharey=True,figsize=(7,6))
    density1 = [pd.concat((density[0][0],density[1][0])),
                pd.concat((density[0][1],density[1][1]))]
    fl = ['(a)','(b)']
    for k11, k1 in enumerate(['N','S']):
        plt.sca(ax[k11])
        density2 = density1[k11]
        if density2.empty:
            continue
        density2['UTbin'] = density2.epochday%1*24//0.5*0.5+0.25
        density2 = density2.groupby('UTbin')['MLT']
        density3 = pd.DataFrame()
        for name, group in density2:
            group1 = group
            if group1.max()-group1.min()>20:
                group1[group1<4] = group1[group1<4]+24
            tmp = pd.DataFrame(
                    {'median':np.median(group1),
                     'p25':np.percentile(group1,25),
                     'p75':np.percentile(group1,75)},
                    index=[name])
            tmp = tmp%24
            density3 = density3.append(tmp)
        plt.plot(density3.index, density3['median'],'ko')
        if k1 == 'S':
            plt.axvline(x=16,color='blue',linestyle='-')
        plt.xlim(0,24)
        plt.xticks(np.arange(0,25,6))
        plt.ylim(0,24)
        plt.yticks(np.arange(0,25,6))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(6))
        plt.grid(dashes=(4,1))
        if k11==1:
            plt.xlabel('UT (hour)')
        plt.ylabel('MLT (hour)')
        plt.text(0.1,0.8,k1,transform=plt.gca().transAxes)
        plt.text(0,1.05,fl[k11], transform=plt.gca().transAxes)
    plt.subplots_adjust(bottom=0.1)
    return


def f3():
    import matplotlib.dates as mdates
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
        rho = [[pd.DataFrame(), pd.DataFrame()],
               [pd.DataFrame(), pd.DataFrame()]]
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
                        print('There is only {:d} '
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

    bae, rho, nsbu= pd.read_pickle(
            os.environ.get('DATAPATH') + 'tmp/w2_f4_02.dat')
    print('Total away and toward days [[AS, AN], [TS, TN]]: \n', nsbu)
    print('Begin figure 1')
    fig, ax = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(6, 8))
    plt.subplots_adjust(
            left=0.13, right=0.82, top=0.93, bottom=0.07,
            wspace=0.11, hspace=0.12)
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
            ll = np.linspace(-3, 3, 11)
            cl = np.arange(-3, 4, 1)
            if k1=='AE':
                ll = np.linspace(0, 300, 11)
                cl = np.arange(0, 301, 100)
            plt.sca(ax[k11, k00])
            baettt = baett.pivot('month', 'uthour', k1)
            # Extend month and ut
            baettt.ix[:, baettt.columns[0]-1] = baettt.ix[:, baettt.columns[-1]]
            baettt.ix[:, baettt.columns[-2]+1] = baettt.ix[:, baettt.columns[0]]
            baettt = baettt.sort_index(axis=1)
            baettt.ix[baettt.index[0]-1, :] = baettt.ix[baettt.index[-1], :]
            baettt.ix[baettt.index[-2]+1, :] = baettt.ix[baettt.index[0], :]
            baettt = baettt.sort_index(axis=0)
            x = baettt.columns # uthour
            y = baettt.index # month
            hc = plt.contourf(x, y, baettt, levels=ll, extend='both',
                              cmap='seismic')
            print('  Average {:s} is: {:5.1f}'.format(k1, baettt.mean().mean()))
            if k1=='Bzm':
                print('  Bz max: {:5.1f}'.format(baettt.max().max()))
                print('  Bz min: {:5.1f}'.format(baettt.min().min()))
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
            plt.tick_params(axis='x', which='major', direction='out', length=4)
            plt.tick_params(axis='y', which='major', direction='out', length=0)
            plt.tick_params(axis='x', which='minor', direction='out', length=2)
            plt.tick_params(axis='y', which='minor', direction='out', length=4)
            plt.text(
                    0.05, 0.8, '('+plb[k11, k00]+')',
                    bbox=dict(facecolor='white', alpha=0.3),
                    transform=plt.gca().transAxes)
            if k11 == 0:
                plt.title(k0)
            if k00 == 0:
                plt.ylabel('Month')
            if k00 == 1:
                axp = plt.gca().get_position()
                cax = plt.axes([0.85, axp.y0, 0.01, axp.height])
                plt.colorbar(hc, cax=cax, orientation='vertical',
                             extend='neither', ticks=cl)
                plt.ylabel(ctt[k11])
                plt.tick_params(axis='both', which='major', length=2)
        for k11, k1 in enumerate(['South', 'North']):
            rhot = rho[k00][k11]
            rhot['month'] = rhot.index.month
            rhot['uthour'] = rhot.index.hour+0.5
            rhott = rhot.groupby(['month', 'uthour']).agg(np.median)
            rhott = rhott.reset_index()
            plt.sca(ax[3+k11, k00])
            rhottt = rhott.pivot('month', 'uthour', 'rrho400')
            # Extend month and ut
            rhottt.ix[:, rhottt.columns[0]-1] = rhottt.ix[:, rhottt.columns[-1]]
            rhottt.ix[:, rhottt.columns[-2]+1] = rhottt.ix[:, rhottt.columns[0]]
            rhottt = rhottt.sort_index(axis=1)
            rhottt.ix[rhottt.index[0]-1, :] = rhottt.ix[rhottt.index[-1], :]
            rhottt.ix[rhottt.index[-2]+1, :] = rhottt.ix[rhottt.index[0], :]
            rhottt = rhottt.sort_index(axis=0)
            hc = plt.contourf(x, y, rhottt, levels=np.linspace(-20, 20, 11),
                              extend='both', cmap='seismic')
            print('  ', k1, ' Density max (%): ', rhottt.max().max())
            print('  ', k1, ' Density min (%): ', rhottt.min().min())
            if k1 is 'South':
                plt.axvline(15+37/60, 0, 1, c='k', ls='--')
            if k1 is 'North':
                plt.axvline(5+25/60, 0, 1, c='k', ls='--')
            plt.xlim(0, 24)
            plt.xticks(np.arange(0, 25, 6))
            plt.ylim(-0.001, 12.001)
            plt.yticks(np.arange(0.5, 12.5),
                       ['', '', 3, '', '', 6, '', '', 9, '', '', 12])
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
            plt.tick_params(axis='x', which='major', direction='out', length=4)
            plt.tick_params(axis='y', which='major', direction='out', length=0)
            plt.tick_params(axis='x', which='minor', direction='out', length=2)
            plt.tick_params(axis='y', which='minor', direction='out', length=4)
            plt.text(
                    0.05, 0.8, '('+plb[k11+3, k00]+')',
                    bbox=dict(facecolor='white', alpha=0.3),
                    transform=plt.gca().transAxes)
            if k00 == 0:
                plt.ylabel('Month')
            if k11 == 1:
                plt.xlabel('UT (hour)')
            if k00 == 1:
                axp = plt.gca().get_position()
                cax = plt.axes([0.85, axp.y0, 0.01, axp.height])
                plt.colorbar(hc, cax=cax, orientation='vertical',
                             extend='neither', ticks=np.arange(-20, 21, 10))
                plt.ylabel(k1[0]+r', $\rho_r$ (%)')
                plt.tick_params(axis='both', which='major', length=2)
    print('End of figure 1\n\n')

    print('Begin figure 2')
    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 4))
    plt.subplots_adjust(
            left=0.13, right=0.95, top=0.90, bottom=0.15,
            wspace=0.11, hspace=0.12)
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
        rhottt.ix[rhottt.index[0]-1, :] = rhottt.ix[rhottt.index[-1], :]
        rhottt.ix[rhottt.index[-2]+1, :] = rhottt.ix[rhottt.index[0], :]
        rhottt = rhottt.sort_index(axis=0)
        hp = plt.plot(rhottt.index, rhottt['median'], 'b')
        print('  Density max (%): ', rhottt['median'].max())
        print('  Density min (%): ', rhottt['median'].min())
        plt.plot(rhottt.index, rhottt.p25,'gray', rhottt.index, rhottt.p75,
                 'gray', linestyle='--')
        plt.axvline(15+37/60, 0, 1, color='k', linestyle='--')
        plt.grid()
        plt.xlim(0, 24)
        plt.xticks(np.arange(0, 25, 6))
        plt.ylim(-30, 30)
        plt.yticks(np.arange(-30, 31, 10))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(6))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.tick_params(axis='both', which='major', length=6)
        plt.tick_params(axis='both', which='minor', length=3)
        if k00 == 0:
            plt.ylabel(r'South, $\rho_r$ (%)')
        plt.xlabel('UT (hour)')
        plt.title(tit)
        plt.text(0.03, 0.91, plb[k00], transform=plt.gca().transAxes)
    print('End of figure 2\n\n')

    print('Begin figure 3')
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
            ((rho.Mlat<=-70) & (rho.Mlat>=-80) &
             (rho.MLT>=11) & (rho.MLT<=13)) ]
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
    print('End of figure 3\n\n')
    return

# END
#-------------------------------------------------------------------------------
if __name__=='__main__':
    plt.close('all')
    a = f3()
    plt.show()
    import gc
    gc.collect()
