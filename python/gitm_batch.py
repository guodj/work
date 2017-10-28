#!/home/guod/anaconda3/bin/python
import gitm
import numpy as np
import matplotlib.pyplot as plt
import gitm_3D_const_alt as g3ca
import os
import pickle
if not os.path.isfile('gitm_batch_variable'):
    gitm_batch_variable={}
    gitm_batch_variable['filename']=''
    gitm_batch_variable['save']='y'
    gitm_batch_variable['savefilename']=''
    gitm_batch_variable['variable']=''
    gitm_batch_variable['log']='n'
    gitm_batch_variable['lineorcontour']='c'
    gitm_batch_variable['constant_surface']='1' # 1, const H; 2, const Lon; 3, const Lat
    gitm_batch_variable['ispolar']='1' # 1, polar; 0, non_polar
    gitm_batch_variable['NS']='1' # 1, north; 0, south
    gitm_batch_variable['min_lat']='50.0' # positive for both poles
    gitm_batch_variable['zmin']='automatic'
    gitm_batch_variable['zmax']='automatic'
    gitm_batch_variable['isvector']='y'
    gitm_batch_variable['neuion']='1' # 1 for neutral wind, 2 for ion drift
    gitm_batch_variable['scale']='-1' # -1 for automatic
    pickle.dump(gitm_batch_variable, open('gitm_batch_variable', 'wb'))

pickle.load(open('gitm_batch_variable', 'rb'))
plt.close('all')
plt.figure(figsize=(9.35, 5.95))
# gitm file name
fn = input('Enter filename to plot[%s]: ' % gitm_batch_variable['filename'])
if len(fn)==0:
    fn =  gitm_batch_variable['filename']
gitm_batch_variable['filename']=fn
# save figure or not
issave = input('Save figure (y or n)? [%s]: ' % gitm_batch_variable['save'])
if len(issave)==0:
    issave = gitm_batch_variable['save']
gitm_batch_variable['save']=issave

g = gitm.GitmBin(fn)
variable_names = (
        'Rho', 'Temperature', 'eTemperature', 'iTemperature',
        'O(!U3!NP)', 'O(!U1!ND)', 'O!D2!N', 'N!D2!N', 'N(!U4!NS)',
        'N(!U2!ND)', 'N(!U2!NP)', 'NO', 'HE', 'H', 'CO!D2!N',
        'V!Dn!N (east)', 'V!Dn!N (north)', 'V!Dn!N (up)',
        'V!Di!N (east)', 'V!Di!N (north)', 'V!Di!N (up)',
        'V!Dn!N (up,O(!U3!NP)           )', 'V!Dn!N (up,O!D2!N              )',
        'V!Dn!N (up,N!D2!N              )', 'V!Dn!N (up,N(!U4!NS)           )',
        'V!Dn!N (up,NO                  )', 'V!Dn!N (up,He                  )',
        'O_4SP_!U+!N', 'O!D2!U+!N', 'N!D2!U+!N', 'N!U+!N', 'NO!U+!N',
        'O(!U2!ND)!U+!N', 'O(!U2!NP)!U+!N', 'H!U+!N', 'He!U+!N', 'e-')
for k00, k0 in enumerate(variable_names):
    print(k00, '     ', k0)
vi = int(input('Enter which variable to plot [%s]: '
                % gitm_batch_variable['variable']))
log10 = input('Enter whether you want log or not (y/n) [n]: ')
log10 = True if log10=='y' else False
lc = input('Enter whether you would like a line plot (l) or contour plot (c) [c]: ')
lc = 'coutour' if lc=='c' else 'line'
vi2 = input('Enter second var to plot as line contour (-1 for None) [-1]: ')
if len(vi2) == 0:
    vi2 = -1
vi2 = int(vi2)
print('1. Constant Altitude Plot')
print('2. Constant Longitude Plot (or Zonal Average)')
print('3. Constant Latitude Plot')
ptype = input('Enter type of plot to make [1]: ')
if len(ptype)==0 or ptype=='1':
    ptype = 'Altitude'
elif ptype=='2':
    ptype = 'Longitude'
elif ptype=='3':
    ptype = 'Latitude'
if ptype == 'Altitude':
    for k0 in range(g.attrs['nAlt']):
        print(k0, '     ', g['Altitude'][0, 0, k0]/1000)
    ialt = int(input('Enter which altitude to plot: '))
    polar = input('Enter polar (1) or non-polar (0) [0]: ')
    if polar=='1':
        polar = 'polar'
        NS=input('Enter North (1) or South (0) [1]: ')
        latmin=input('Enter minimum latitude to plot [50.0] :')
        if len(latmin)==0:
            latmin=50
        latmin = abs(float(latmin))
        if NS=='1' or len(NS)==0:
            nlat, slat = 90, latmin
        elif NS=='0':
            nlat, slat = -latmin, -90
        useLT=True
        dlonlt=6
        dlat=10
        ax = plt.subplot(111, polar=True)
    elif polar=='0' or len(polar)==0:
        polar = 'rectangular'
        nlat, slat = 90, -90
        ax = plt.subplot(111)
        dlonlt=30
        dlat=30
        useLT=False
    zmin = input('Enter minimum [automatic]: ')
    if len(zmin)==0:
        zmin=None
    else:
        zmin=float(zmin)
    zmax = input('Enter maximum [automatic]: ')
    if len(zmax)==0:
        zmax=None
    else:
        zmax=float(zmax)
    axt, hc = g3ca.contour_single(
            ax, vl[vi], polar, g, ialt=ialt, nlat=nlat, slat=slat,
            dlonlt=dlonlt, zmax=zmax, zmin=zmin, log10=log10,
            useLT=useLT, dlat=dlat)
    plt.colorbar(hc, ax=axt, pad=0.1)
    if vi2!=-1:
        axt, hc = g3ca.contour_single(
                ax, vl[vi2], polar, g, ialt=ialt, nlat=nlat, slat=slat,
                dlonlt=dlonlt, zmax=zmax, zmin=zmin, fill=False, colors='k',
                alpha=0.5, dlat=dlat, useLT=useLT)
        plt.clabel(hc, fmt='%.1f')
    vectoryn = input('Enter whether you want vectors or not (y/n) [y]: ')
    if vectoryn=='y' or len(vectoryn)==0:
        neuion = input('Enter plot neutral winds (1) or ions (0) [1]:')
        if neuion=='1' or len(neuion)==0:
            neuion = 'neutral'
        elif neuion=='0':
            neuion = 'ion'
        scalel=[None, 10, 50, 100, 200, 250, 500, 750, 1000, 1500, 2000, 3000]
        print('-1   :   automatic selection')
        for k00, k0 in enumerate(scalel[1:]):
            print(k00+1, '    :   ', scalel[k00+1])
        scalew = input('Enter velocity factor [-1]:')
        if scalew=='-1' or len(scalew)==0:
            scalew = 0
        scalew = int(scalew)
        axt, hv = g3ca.vector_single(
                ax, g, neuion, polar, ialt=ialt, nlat=nlat, slat=slat,
                scale=scalel[scalew], alpha=0.6, useLT=useLT,
                dlat=dlat, dlonlt=dlonlt)
    plt.title(g[vl[vi]].attrs['name']+' at '+
              '%.1f'%(g['Altitude'][0, 0, ialt]/1000)+
              ' km Altitude at '+ g['time'].strftime('%d-%b-%y %H:%M')+ ' UT')
if ptype == 'Longitude':
    print('Unimplemented')
if ptype == 'Latitude':
    print('Unimplemented')
plt.savefig(fnpdf)
