#--------------------------------------------------------------------------------
# Functions to handle Omni data.
#
# By Dongjie, UM, on Thu Sep 22 01:57:23 CST 2016
#
# Contain:
#     get_omni: Get Omni data.
#     plot_omni: Plot Omni data
#--------------------------------------------------------------------------------

# Global imports
import pandas as pd
import numpy as np
import os

OMNIDATADIR = '/home/guod/data/omni/'
def get_omni_backup(bdate, edate,variables,res='1h'):
    global OMNIDATADIR
    # Get multiple days omni data.
    # variables should be list or tuple even if only one variable.
    # This function is replaced by get_omni() for better time performance
    bdate = pd.Timestamp(bdate)
    edate = pd.Timestamp(edate)
    years = np.arange(bdate.year, edate.year+1)
    variables = list(variables) # I want to apppend date and time.
    def parser_low(date):
            return pd.to_datetime(date,format='%Y %j %H')

    def parser_high(date):
            return pd.to_datetime(date,format='%Y %j %H %M')
    # low resolution
    if res in ['1day', '1d', '1D','1hour', '1Hour', '1h', '1H']:
        fp = 'low'
        column_name = [
                'year','day','hour',
                't1','t2','t3','t4','t5','t6','t7','t8','t9',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't10','t11','t12','t13','t14',
                'ProTmp','ProDen', 'V',
                't15','t16','t17','t18','t19','t20','t21','t22','t23','t24','t25','t26','t27',
                'Kp','R','DST','AE',
                't28','t29','t30','t31','t32','t33','t34',
                'ap','f107','PC','AL','AU',
                't35']
        na_value={'Bx':[999.9],'Bye':[999.9],'Bze':[999.9],'Bym':[999.9],'Bzm':[999.9],
                   'ProTmp':[9999999.],'ProDen':[999.9],'V':[9999.],
                   'Kp':[99],'R':[999],'DST':[99999],'AE':[9999]}
        if res in ['1day', '1d', '1D']:
            # Use list in order to keep consistency with other resolutions
            fname = [OMNIDATADIR+'low_res_omni/omni_01_av.dat']
        if res in ['1hour', '1Hour', '1h', '1H']:
            fname = [OMNIDATADIR+'low_res_omni/omni2_{:4d}.dat'.format(k) for k in years]
        index_clm = ['year','day','hour']
    # high resolution
    if res in ['5minute', '5minutes', '5m','1minute', '1m']:
        fp = 'high'
        column_name = [
                'year','day','hour','minute',
                't1','t2','t3','t4','t5','t6','t7','t8','t9','t10',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't11','t12',
                'V','Vx','Vy','Vz','ProDen','ProTmp',
                't13','t14','t15','t16','t17','t18','t19','t20','t21','t22',
                'AE','AL','AU','SYMD','SYMH','ASYD','ASYH','PC',
                't23']
        na_value={'Bx':[9999.99],'Bye':[9999.99],'Bze':[9999.99],'Bym':[9999.99],'Bzm':[9999.99],
                   'V':[99999.9],'Vx':[99999.9],'Vy':[99999.9],'Vz':[99999.9],
                   'ProTmp':[9999999.],'ProDen':[999.99],
                   'AE':[99999],'AL':[99999],'AU':[99999],
                   'SYMD':[99999],'SYMH':[99999],'ASYD':[99999],'ASYH':[99999],'PC':[999.99]}
        if res in ['5minute', '5minutes', '5m']:
            fname =[OMNIDATADIR+'high_res_omni/omni_5min{:4d}.asc'.format(k) for k in years]
            column_name.extend(['t24','t25','t26'])
        if res in ['1minute', '1m']:
            fname =[OMNIDATADIR+'high_res_omni/omni_min{:4d}.asc'.format(k) for k in years]
        index_clm = ['year','day','hour','minute']
    variables.extend(index_clm)
    parser = parser_low if fp is 'low' else parser_high
    omni_data = [
            pd.read_csv(
                    fn, delim_whitespace=True,
                    header=None, names=column_name, usecols=variables,
                    index_col='date', squeeze=True,
                    parse_dates={'date':index_clm}, date_parser=parser,
                    na_values=na_value)
            for fn in fname if os.path.isfile(fn)]
    if omni_data:
        omni_data = pd.concat(omni_data)
        omni_data = omni_data[bdate:edate]
        return omni_data
    else:
        return pd.DataFrame()

def get_omni(bdate, edate,variables,res='1m'):
    # Get multiple days omni data.
    # variables should be list or tuple even if only one variable.
    import datetime as dt
    global OMNIDATADIR
    bdate = pd.Timestamp(bdate)
    edate = pd.Timestamp(edate)
    years = np.arange(bdate.year, edate.year+1)
    # low resolution
    if res in ['1day', '1d', '1D','1hour', '1Hour', '1h', '1H']:
        fp = 'low'
        column_name = [
                'year','day','hour',
                't1','t2','t3','t4','t5','t6','t7','t8','t9',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't10','t11','t12','t13','t14',
                'ProTmp','ProDen', 'V',
                't15','t16','t17','t18','t19','t20','t21','t22','t23','t24','t25','t26','t27',
                'Kp','R','DST','AE',
                't28','t29','t30','t31','t32','t33','t34',
                'ap','f107','PC','AL','AU',
                't35']
        na_value={'Bx':[999.9],'Bye':[999.9],'Bze':[999.9],'Bym':[999.9],'Bzm':[999.9],
                   'ProTmp':[9999999.],'ProDen':[999.9],'V':[9999.],
                   'Kp':[99],'R':[999],'DST':[99999],'AE':[9999]}
        if res in ['1day', '1d', '1D']:
            # Use list in order to keep consistency with other resolutions
            fname = [OMNIDATADIR+'low_res_omni/omni_01_av.dat']
            freq = 'D'
        if res in ['1hour', '1Hour', '1h', '1H']:
            fname = [OMNIDATADIR+'low_res_omni/omni2_{:4d}.dat'.format(k) for k in years]
            freq = 'H'
    # high resolution
    if res in ['5minute', '5minutes', '5m','1minute', '1m']:
        fp = 'high'
        column_name = [
                'year','day','hour','minute',
                't1','t2','t3','t4','t5','t6','t7','t8','t9','t10',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't11','t12',
                'V','Vx','Vy','Vz','ProDen','ProTmp',
                't13','t14','t15','t16','t17','t18','t19','t20','t21','t22',
                'AE','AL','AU','SYMD','SYMH','ASYD','ASYH','PC',
                't23']
        na_value={'Bx':[9999.99],'Bye':[9999.99],'Bze':[9999.99],'Bym':[9999.99],'Bzm':[9999.99],
                   'V':[99999.9],'Vx':[99999.9],'Vy':[99999.9],'Vz':[99999.9],
                   'ProTmp':[9999999.],'ProDen':[999.99],
                   'AE':[99999],'AL':[99999],'AU':[99999],
                   'SYMD':[99999],'SYMH':[99999],'ASYD':[99999],'ASYH':[99999],'PC':[999.99]}
        if res in ['5minute', '5minutes', '5m']:
            fname =[OMNIDATADIR+'high_res_omni/omni_5min{:4d}.asc'.format(k) for k in years]
            column_name.extend(['t24','t25','t26'])
            freq = '5T'
        if res in ['1minute', '1m']:
            fname =[OMNIDATADIR+'high_res_omni/omni_min{:4d}.asc'.format(k) for k in years]
            freq = 'T'
    omni_data = []
    for k0, k1 in zip(years, fname):
        if not os.path.isfile(k1):
            continue
        tmp = pd.read_csv(
                k1, delim_whitespace=True,
                header=None, names=column_name, usecols=variables,
                squeeze=True,
                na_values=na_value)
        if not tmp.empty:
            tmp = pd.DataFrame(tmp)
            tmp['date'] = pd.date_range(
                    dt.datetime(k0, 1, 1, 0, 0, 0),
                    dt.datetime(k0, 12, 31, 23, 59, 0), freq=freq)
            tmp.set_index('date', drop=True, inplace=True)
            omni_data.append(tmp)
    if omni_data:
        omni_data = pd.concat(omni_data)
        omni_data = omni_data[bdate:edate]
        return omni_data # datatype: a pd.DataFrame, not a pd.Series
    else:
        return pd.DataFrame()

#----------------------------------------
def plot_omni(ax, bdate, edate, variables, res='1hour', **kwargs):
    # ax should be a list of axes
    # variables should be a list or tuple even if only one variable
    # kwargs are parameters for plot
    # Plot variables in panels
    # res can be '1D','1H','5m','1m'
    #----------------------------------------
    # variablename gives the column names in pd.Dataname and the possible names
    # used in ticklables
    import matplotlib.pyplot as plt
    variablename = {'Bx':'$B_x$ (nT)',
                    'Bye':'GSE $B_y$ (nT)', 'Bym':'GSM $B_y$ (nT)',
                    'Bze':'GSE $B_z$ (nT)', 'Bzm':'GSM $B_z$ (nT)',
                    'AE':'AE (nT)', 'AL':'AL (nT)', 'AU':'AU (nT)',
                    'PC':'PC',
                    'V':'$V_SW$ (km/s)', 'Vx':'$V_x$ (km/s)',
                    'Vy':'$V_y$ (km/s)', 'Vz':'$V_z$ (km/s)', # flow speed
                    'ProDen':'Proton Density ($N/cm^3$)',
                    'ProTmp':'Proton Temperature (k)',
                    'SYMD': 'SYM/D (nT)','SYMH': 'SYM/H (nT)',
                    'ASYD': 'ASY/D (nT)','ASYH': 'ASY/H (nT)',
                    'Kp':'Kp', 'R':'R','DST':'DST','ap':'ap','f107':'f10.7'}
    omni_data = get_omni(bdate,edate,variables,res)
    omni_data = pd.DataFrame(omni_data)
    for k00, k0 in enumerate(variables):
        plt.sca(ax[k00]) if len(variables)>1 else plt.sca(ax)
        plt.plot(omni_data.index, omni_data[k0], **kwargs)
        plt.ylabel(variablename[k0])
    return
#END
#--------------------------------------------------------------------------------
# TEST
if __name__ == '__main__':
    # Test plot_omni
    #    fig,ax = plt.subplots(1,1,sharex=True)
    #    plot_omni(ax,'2010-1-1','2010-1-31',['Bx'],res='5minute')
    # Test get_omni and get_omni2
    from timeit import time
    import matplotlib.pyplot as plt
    begintime = time.time()
    omni1 = get_omni_backup('2010-1-1', '2010-12-31', ['Bx'], res='1m')
    endtime = time.time()
    print(endtime-begintime)
    begintime = time.time()
    omni2 = get_omni('2010-1-1', '2010-12-31', ['Bx'], res='1m')
    endtime = time.time()
    print(endtime-begintime)
    plt.plot(omni1-omni2['Bx'])

