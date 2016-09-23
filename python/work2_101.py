#--------------------------------------------------------------------------------
# For the second project
#
# By Dongjie, USTC, Sat Sep 17 09:32:19 CST 2016
#
#--------------------------------------------------------------------------------

#Global imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pdb   # set breakpoint
from champ_grace import *
from omni import *



def get_sblist():
    """ Get solar wind sector polarity reversing dates

    Args:
        No input
    Returns:
        A dataframe indexed by dates with columns sbtype, lday, rday and season.
    """
    sb_fname = '/data/SBlist/SBlist.txt'
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


def f1():
    # imf and indices variation during epoch days of sblist
    sblist = get_sblist()
    sblist = sblist['2002-1-1':'2010-12-31']
    data = [pd.DataFrame(),pd.DataFrame()]
    if False:  # data preparation
        for k00, k0 in enumerate(['away-toward', 'toward-away']):
            sbtmp = sblist[sblist.sbtype==k0]
            for k1 in sbtmp.index:
                print(k1,k0)
                bdate = k1-pd.Timedelta('3D')
                edate = k1+pd.Timedelta('2D')
                data_tmp = get_omni(bdate,edate,['Bx','Bym','Bzm','AE'],res='1h')
                #print(data_tmp)
                data_tmp['epochday'] = (data_tmp.index - k1)/pd.Timedelta('1D')
                data[k00] = data[k00].append(data_tmp)
        pd.to_pickle(data,'/data/tmp/w2_07.dat')
    # END of data preperation

    data = pd.read_pickle('/data/tmp/w2_07.dat')
    datagroup = [
            data[k].groupby([data[k].index.month,np.floor(data[k].epochday*24)])
            for k in [0,1]]
    datagroup = [datagroup[k].median() for k in [0,1]]
    for k in [0,1]:
        datagroup[k].index.names = ('month', 'epochhour')
        datagroup[k] = datagroup[k].reset_index().pivot(index='epochhour', columns='month')
    fig,ax = plt.subplots(4,2,sharex=True,sharey=True,figsize=(7.76,8))
    for k in range(2):
        for k11,k1 in enumerate(['Bx','Bym','Bzm','AE']):
            tl = ['Bx (nT)','By (nT)','Bz (nT)','AE']
            levels = np.linspace(0,400,11) if k1 is 'AE' else np.linspace(-4,4,11)
            ticks = np.arange(0,401,100) if k1 is 'AE' else np.arange(-4,5,2)
            plt.sca(ax[k11,k])
            data1 = datagroup[k][k1]
            hc2 = plt.contourf(data1.columns, data1.index/24, data1.values,levels=levels,
                    cmap='bwr')
            plt.xlim([1,12])
            plt.xticks(np.arange(1,13))
            plt.ylim([-3,3])
            plt.yticks(np.arange(-3,4,1),fontsize=14)
            plt.tick_params(axis='both',which='major',direction='out',length=5)
            plt.tick_params(axis='both',which='minor',direction='out',length=3)
            if k is 1:
                axpo = np.array(plt.gca().get_position())
                cax = plt.gcf().add_axes((axpo[1,0]+0.005,axpo[0,1],0.01,axpo[1,1]-axpo[0,1]))
                cbar = plt.colorbar(mappable=hc2,cax=cax,ticks=ticks)
                cbar.set_label(tl[k11])
                plt.tick_params('both',length=4)
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
    so it is a good location for the research of UT variation
    """
    def percentile(n):
        def percentile_(x):
            return np.percentile(x,n)
        percentile_.__name__ = 'percentile_%s' % n
        return percentile_

    sblist = get_sblist()
    sblist = sblist['2001-1-1':'2010-12-31']
    density = [[pd.DataFrame(),pd.DataFrame()],[pd.DataFrame(),pd.DataFrame()]] # [[ATN,ATS],[TAN,TAS]]
    sbn = [0,0]
    if False:
        for k00,k in enumerate(['away-toward','toward-away']):
            sbtmp = sblist[sblist.sbtype==k]
            for k1 in sbtmp.index:
                #for k2 in ['champ','grace']:
                for k2 in ['grace']:  # only consider the grace
                    rho = get_champ_grace_data(k1+pd.TimedeltaIndex(range(-3,3),'D'), k2)
                    if rho.empty:
                        print('no data around',k1)
                        continue
                    # The geomagnetic activity is not considered, because cusp density
                    # enhancement occurs in both geomagnetic quiet and active conditions
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
                    rho1 = rho[rho.lat3==90]  # only consider the grace
                    #Make sure that there is data in each day
                    if len(np.unique(rho1.index.dayofyear))!=6:
                        print('there is data gap around', k1)
                        continue
                    rho1['epochday'] = (rho1.index-k1)/pd.Timedelta('1D')
                    rho1['rrho400'] = 100*(rho1.rho400-rho1['rho400'].mean())/rho1['rho400'].mean()
                    rho1['rrho'] = 100*(rho1.rho-rho1['rho'].mean())/rho1['rho'].mean()
                    density[k00][0] = density[k00][0].append(
                            rho1[['epochday','MLT','Mlat','rho','rrho','rho400','rrho400']])

                    #rho2 = rho[rho.lat3<=-87]  # south pole
                    rho2 = rho[rho.lat3==-90]  # only consider the grace
                    if len(np.unique(rho2.index.dayofyear))!=6:
                        print('there is data gap around', k1)
                        continue
                    rho2['epochday'] = (rho2.index-k1)/pd.Timedelta('1D')
                    rho2['rrho400'] = 100*(rho2.rho400-rho2['rho400'].mean())/rho2['rho400'].mean()
                    rho2['rrho'] = 100*(rho2.rho-rho2['rho'].mean())/rho2['rho'].mean()
                    density[k00][1] = density[k00][1].append(
                            rho2[['epochday','MLT','Mlat','rho','rrho','rho400','rrho400']])
                    sbn[k00] = sbn[k00]+1
                    print(sbn)
        pd.to_pickle(density, '/data/tmp/w2_13.dat')
    # END of data preperation

    # Pole density variation as a function of epoch time at different seasons and sbtype.
    density = pd.read_pickle('/data/tmp/w2_13.dat')
    fig,ax = plt.subplots(4,4,sharex=True,sharey=True,figsize=(8,8))
    fl = [['(a1)','(a2)','(a3)','(a4)'],['(b1)','(b2)','(b3)','(b4)'],
          ['(c1)','(c2)','(c3)','(c4)'],['(d1)','(d2)','(d3)','(d4)']]
    # case number in each season catagary
    nn = np.zeros([4,4])*np.nan
    for k00,k in enumerate(['away-toward','toward-away']):
        for k11, k1 in enumerate(['N','S']):
            density1 = density[k00][k11]
            #density1 = density1['2002-1-1':'2004-12-31']
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
                density2 = density1[fp]
                nn[k22, k00*2+k11] = len(np.unique(
                        density2[(density2.epochday>=0) & (density2.epochday<1)].index.date))
                density2['epochbin'] = density2.epochday*24//1.5*1.5+0.75
                density2 = density2.groupby('epochbin')['rrho400'].agg(
                        [np.median, percentile(25),percentile(75)])
                density2.columns = ['median', 'p25', 'p75']
                plt.plot(density2.index/24, density2['p25'],'gray',
                         density2.index/24, density2['p75'],'gray',
                         linestyle='--',dashes=(2,1),linewidth=1)
                plt.plot(density2.index/24, density2['median'],'b',linewidth=2)
                plt.xlim(-3,3)
                plt.xticks(np.arange(-3,4,1),('',-2,-1,0,1,2,''))
                #plt.gca().xaxis.set_minor_locator(AutoMinorLocator(4))
                if k1 is 'S':
                    plt.vlines(np.arange(-3,3)+15.5/24,-30,60,linestyle='--',linewidth=1,color='r')
                if k1 is 'N':
                    plt.vlines(np.arange(-3,3)+5.5/24,-30,60,linestyle='--',linewidth=1,color='r')
                plt.ylim(-30,60)
                plt.yticks(np.arange(-30,61,30))
                #plt.grid(which='minor',dashes=(4,1))
                #plt.grid(which='major',axis='y',dashes=(4,1))
                plt.grid(axis='y',dashes=(4,1))
                plt.tick_params(axis='both',which='major',direction='out',length=5)
                if k00*2+k11==0:
                    plt.ylabel(r'$\Delta\rho$ (%)')
                if k22==3:
                    plt.xlabel('Epoch Time (day)',fontsize=12)
                plt.text(0.1,0.8,k1,transform=plt.gca().transAxes)
                #plt.text(0,1.05,fl[k22][k00*2+k11], transform=plt.gca().transAxes)
    plt.subplots_adjust(left=0.1,wspace=0.04)
    plt.text(0.21,0.94,'Away - Toward',transform=plt.gcf().transFigure)
    plt.text(0.61,0.94,'Toward - Away',transform=plt.gcf().transFigure)
    plt.text(0.91,0.8,'Feb - Apr',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.59,'Aug - Oct',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.38,'May - Jul',transform=plt.gcf().transFigure,fontsize=11)
    plt.text(0.91,0.17,'Nov - Jan',transform=plt.gcf().transFigure,fontsize=11)
    print(nn)

    # Density variations at solar maximum and minimum.
    fig,ax = plt.subplots(2,1,sharex=True,sharey=True,figsize=(6.35,7.02))
    density1 = density[0][1] # for away-toward and south pole
    fl = ['(a)','(b)']
    for k00,k0 in enumerate(['Solar maximum','Solar minimum']):
        plt.sca(ax[k00])
        if k0 is 'Solar maximum':
            density2 = density1['2002-1-1':'2004-12-31']
        if k0 is 'Solar minimum':
            density2 = density1['2008-1-1':'2010-12-31']
        density2 = density2[(density2.index.month>=8) & (density2.index.month<=10)]
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
        plt.vlines(np.arange(-3,3)+15.5/24,-30,61,linestyle='--',linewidth=1,color='r')
        plt.ylim(-30,60)
        plt.yticks(np.arange(-30,61,30))
        plt.tick_params(axis='both',which='major',direction='out',length=5)
        #plt.grid(which='minor',dashes=(4,1))
        #plt.grid(which='major',axis='y',dashes=(4,1))
        plt.grid(axis='y',dashes=(4,1))
        plt.ylabel(r'$\Delta\rho$ (%)',fontsize=14)
        if k00==1:
            plt.xlabel('Epoch Time (day)',fontsize=14)
        if k00==0:
            a = plt.title('Year: 02 - 04')
        if k00==1:
            a = plt.title('Year: 08 - 10')
        a.set_position((0.5,1.06))
        plt.text(0.1,0.8,'S',transform=plt.gca().transAxes)
        plt.text(0,1.05,fl[k00], transform=plt.gca().transAxes)
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
# END
#--------------------------------------------------------------------------------
if __name__=='__main__':
    plt.close('all')
    a = f1()
    plt.show()
    import gc
    gc.collect()
