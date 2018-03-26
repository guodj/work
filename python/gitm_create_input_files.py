#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie, 2016-12-08 10:04, UM
# Comment: create input files for Gitm
#-------------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------
#                              Case 1
#------------------------------------------------------------
# Basic information
homepath = os.environ.get('HOME')
btime = pd.Timestamp('2003-03-19 00:00:00')
etime = pd.Timestamp('2003-03-24 00:00:00')
dt = pd.Timedelta('1min')
time = pd.date_range(btime, etime, freq=dt)
#----------------------------------------
# Write imf.dat
dim = time.shape
imfby1 = np.ones(dim)*0 # run1
reversalt = etime-pd.Timedelta('6hour')
mm = time>reversalt
imfby2 = imfby1.copy()
imfby2[mm] = imfby2[mm]+4.9
imfbz1 = np.ones(dim)*-5
imfbz2 = imfbz1.copy()
imfbz2[mm] = imfbz2[mm]-1.1
imfbx = np.zeros(dim)
vx = np.ones(dim)*-450
vy = np.zeros(dim)
vz = np.zeros(dim)
den = np.ones(dim)*5
tem = np.ones(dim)*50000
fimf1 = open(homepath+'/tmp/imf1.dat', 'w')
fimf1.write('#START\n')
for k00, k0 in enumerate(time):
    # Date and Time
    fimf1.write('{:5d}{:3d}{:3d}{:3d}{:3d}{:3d}{:4d}'.format(
            k0.year, k0.month, k0.day, k0.hour,
            k0.minute, k0.second, k0.microsecond))
    # IMF Bx, By, Bz
    fimf1.write('{:8.2f}{:8.2f}{:8.2f}'.format(
            imfbx[k00], imfby1[k00], imfbz1[k00]))
    # Solar wind velocity, Vx, Vy, Vz
    fimf1.write('{:9.2f}{:9.2f}{:9.2f}'.format(
            vx[k00], vy[k00], vz[k00]))
    # Solar wind density and temperature
    fimf1.write('{:9.2f}{:13.1f}'.format(den[k00], tem[k00]))
    fimf1.write('\n')
fimf1.close()

fimf2 = open(homepath+'/tmp/imf2.dat', 'w')
fimf2.write('#START\n')
for k00, k0 in enumerate(time):
    # Date and Time
    fimf2.write('{:5d}{:3d}{:3d}{:3d}{:3d}{:3d}{:4d}'.format(
            k0.year, k0.month, k0.day, k0.hour,
            k0.minute, k0.second, k0.microsecond))
    # IMF Bx, By, Bz
    fimf2.write('{:8.2f}{:8.2f}{:8.2f}'.format(
            imfbx[k00], imfby2[k00], imfbz2[k00]))
    # Solar wind velocity, Vx, Vy, Vz
    fimf2.write('{:9.2f}{:9.2f}{:9.2f}'.format(
            vx[k00], vy[k00], vz[k00]))
    # Solar wind density and temperature
    fimf2.write('{:9.2f}{:13.1f}'.format(den[k00], tem[k00]))
    fimf2.write('\n')
fimf2.close()
#----------------------------------------
# Write sme.dat
au = np.ones(dim)*50+np.round(np.random.normal(0, 5, dim))
al = np.ones(dim)*-50+np.round(np.random.normal(0, 5, dim))
ae = au-al
t1 = np.zeros(dim)
t2 = np.zeros(dim)
t3 = np.zeros(dim)
t4 = np.zeros(dim)
t5 = np.zeros(dim)
fsme = open(homepath+'/tmp/sme.dat', 'w')
for k00, k0 in enumerate(time):
    # Date and time
    fsme.write('{:5d}{:3d}{:3d}{:3d}{:3d}{:3d}'.format(
            k0.year, k0.month, k0.day, k0.hour, k0.minute, k0.second))
    # AE, AL, AU
    fsme.write('{:6.0f}{:6.0f}{:6.0f}'.format(
            ae[k00], al[k00], au[k00]))
    # Something not important
    fsme.write('{:6.0f}{:6.0f}{:6.0f}{:6.0f}{:6.0f}'.format(
            t1[k00], t2[k00], t3[k00], t4[k00], t5[k00]))
    fsme.write('\n')
fsme.close()
#----------------------------------------
# Write onsets.dat
time1 = btime-pd.Timedelta('24hour')
time2 = etime+pd.Timedelta('24hour')
timetime =[time1, time2]
fonsets = open(homepath+'/tmp/onsets.dat', 'w')
for k in range(2):
    timek = timetime[k]
    # Date and time
    fonsets.write('{:5d}{:3d}{:3d}{:3d}{:3d}'.format(
        timek.year, timek.month, timek.day, timek.hour, timek.minute))
    # Onset latitude and local time
    fonsets.write('{:6.2f}{:6.2f}'.format(70, 22))
    fonsets.write('\n')
fonsets.close()
#----------------------------------------
# display results
# Run1
fig1, ax1 = plt.subplots(3, 1, sharex=True, figsize=(7, 8))
plt.sca(ax1[0])
plt.title('Run1')
plt.plot(time, imfby1, 'b', time, imfbz1, 'r')
hline = [5, -5, -1]
[plt.axhline(k, linestyle='--', color='gray') for k in hline]
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel(r'IMF $B_y$, $B_z$')

plt.sca(ax1[1])
plt.plot(time, vx, 'k')
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel(r'$V_x$')

plt.sca(ax1[2])
plt.plot(time, au, 'gray', time, al, 'gray', time, ae, 'k')
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel('AE, AU, AL')

xticks = pd.date_range(btime, etime, freq='6h')
xticklabels = ((xticks-btime).total_seconds()/3600).astype(np.int)
plt.xticks(xticks, xticklabels)
plt.xlabel('Hours from '+btime.strftime('%Y-%m-%d %H:%M:%S'))
# Run2
fig2, ax2 = plt.subplots(3, 1, sharex=True, figsize=(7, 8))
plt.sca(ax2[0])
plt.title('Run2')
plt.plot(time, imfby2, 'b', time, imfbz2, 'r')
hline = [5, -5, -1]
[plt.axhline(k, linestyle='--', color='gray') for k in hline]
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel(r'IMF $B_y$, $B_z$')

plt.sca(ax2[1])
plt.plot(time, vx, 'k')
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel(r'$V_x$')

plt.sca(ax2[2])
plt.plot(time, au, 'gray', time, al, 'gray', time, ae, 'k')
vline = pd.date_range(btime, etime, freq='1D')
[plt.axvline(k, linestyle='--', color='gray') for k in vline]
plt.ylabel('AE, AU, AL')

xticks = pd.date_range(btime, etime, freq='6h')
xticklabels = ((xticks-btime).total_seconds()/3600).astype(np.int)
plt.xticks(xticks, xticklabels)
plt.xlabel('Hours from '+btime.strftime('%Y-%m-%d %H:%M:%S'))
plt.show()
# END OF CASE 1
