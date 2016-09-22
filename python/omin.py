#--------------------------------------------------------------------------------
# Functions to handle Omin data.
#
# By Dongjie, USTC, on Thu Sep 22 01:57:23 CST 2016
#
# Contain:
#     get_omin: Get Omin data.
#--------------------------------------------------------------------------------

# Global imports
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# VariableName gives the column names in pd.Dataname and the possible names
# used in ticklables
VARIABLENAME = {'Bx':'$B_x$ (nT)',
                'Bye':'GSE $B_y$ (nT)', 'Bym':'GSM $B_y$ (nT)',
                'Bze':'GSE $B_z$ (nT)', 'Bzm':'GSM $B_z$ (nT)',
                'AE':'AE (nT)', 'AL':'AL (nT)', 'AU':'AU (nT)',
                'PC':'PC',
                'V':'$V_SW$ (km/s)', 'Vx':'$V_x$ (km/s)', 'Vy':'$V_y$ (km/s)', 'Vz':'$V_z$ (km/s)', # flow speed
                'ProDen':'Proton Density ($N/cm^3$)', 'ProTmp':'Temperature (k)',
                'SYMD': 'SYM/D (nT)','SYMH': 'SYM/H (nT)','ASYD': 'ASY/D (nT)','ASYH': 'ASY/H (nT)',
                'Kp':'Kp', 'R':'R','DST':'DST','ap':'ap','f107':'f10.7'}

def get_omin(bdate, edate,variables,res='1h'):
    # variables should be list or tuple even if only one variable.
    bdate = pd.Timestamp(bdate)
    edate = pd.Timestamp(edate)
    years = np.arange(bdate.year, edate.year+1)
    # low resolution
    if res in ['1day', '1d', '1D','1hour', '1Hour', '1h', '1H']:
        fp = 'low'
        column_name = (
                'year','day','hour',
                't1','t2','t3','t4','t5','t6','t7','t8','t9',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't10','t11','t12','t13','t14',
                'ProTmp','ProDen', 'V',
                't15','t16','t17','t18','t19','t20','t21','t22','t23','t24','t25','t26','t27',
                'Kp','R','DST','AE',
                't28','t29','t30','t31','t32','t33','t34',
                'ap','f107','PC','AL','AU',
                't35')
        na_values={'Bx':[999.9],'Bye':[999.9],'Bze':[999.9],'Bym':[999.9],'Bzm':[999.9],
                   'ProTmp':[9999999.],'ProDen':[999.9],'V':[9999.],
                   'Kp':[99],'R':[999],'DST':[99999],'AE':[9999]}
        if res in ['1day', '1d', '1D']:
            # Use list in order to keep consistency with other resolutions
            fname = ['/data/omni/low_res_omni/omni_01_av.dat']
        if res in ['1hour', '1Hour', '1h', '1H']:
            fname = ['/data/omni/low_res_omni/omni2_{:4d}.dat'.format(k) for k in years]
    # high resolution
    if res in ['5minute', '5minutes', '5m''1minute', '1m']:
        fp = 'high'
        column_name = (
                'year','day','hour','minute',
                't1','t2','t3','t4','t5','t6','t7','t8','t9','t10',
                'Bx', 'Bye','Bze','Bym','Bzm',
                't11','t12',
                'V','Vx','Vy','Vz','ProDen','ProTmp',
                't13','t14','t15','t16','t17','t18','t19','t20','t21','t22',
                'AE','AL','AU','SYMD','SYMH','ASYD','ASYH','PC',
                't23')
        na_values={'Bx':[9999.99],'Bye':[9999.99],'Bze':[9999.99],'Bym':[9999.99],'Bzm':[9999.99],
                   'V':[99999.9],'Vx':[99999.9],'Vy':[99999.9],'Vz':[99999.9],
                   'ProTmp':[9999999.],'ProDen':[999.99],
                   'AE':[99999],'AL':[99999],'AU':[99999],
                   'SYMD':[99999],'SYMH':[99999],'ASYD':[99999],'ASYH':[99999],'PC':[999.99]}
        if res in ['5minute', '5minutes', '5m']:
            fname =['/data/omni/high_res_omni/omni_5min{:4d}.asc'.format(k) for k in years]
        if res in ['1minute', '1m']:
            fname =['/data/omni/high_res_omni/omni_min{:4d}.asc'.format(k) for k in years]
    omni_data = [
            pd.read_csv(
                    fn, delim_whitespace=True,
                    header=None, names=column_name, index_col=False, usecols=variables,squeeze=True,
                    na_values=na_values)
            for fn in fname if os.path.isfile(fn)]
    if omni_data:
        omni_data = pd.concat(omni_data)
        if fp is 'low':
            omni_data['mintute']=0
        omni_data['date'] = pd.to_datetime(omni_data[['year','day','hour']])
    return fname
#END
#--------------------------------------------------------------------------------
# TEST
if __name__ == '__main__':
    a = get_omin('2005-1-1','2005-1-10',['ProTmp'],res='1h')

