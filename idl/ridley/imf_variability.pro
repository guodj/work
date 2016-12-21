filename = 'imf_by-0+0.0_bz-5+0.0_vx-400+100.dat'

bx_mean = 0.0 
bx_amp  = 0.0

by_mean =  0.0
by_amp  =  0.0

bz_mean =   -5.0
bz_amp  =    0.0

vx_mean = -400.0
vx_amp  =   100.0

vy_mean =   0.0
vy_amp  =   0.0

vz_mean =   0.0
vz_amp  =   0.0

n_mean =  5.0
n_amp  =  0.0

itime = [2011, 3, 19, 0, 0, 0]
c_a_to_r, itime, basetime

dt = 1.0 * 60.0

nDays  = 4
nTimes = nDays * (24.0*3600.0)/dt + 1

time = findgen(nTimes) * dt / 60.0

bx = bx_mean + bx_amp * randomn(seed,nTimes)
by = by_mean + by_amp * randomn(seed,nTimes)
bz = bz_mean + bz_amp * randomn(seed,nTimes)
vx = vx_mean + vx_amp * randomn(seed,nTimes)
vy = vy_mean + vy_amp * randomn(seed,nTimes)
vz = vz_mean + vz_amp * randomn(seed,nTimes)
n  =  n_mean +  n_amp * randomn(seed,nTimes)

close,/all
openw,1,filename

printf,1,''
printf,1,'#START'
for i=0,nTimes-1 do begin
    rtime = basetime + time(i)*60.0
    c_r_to_a, itime, rtime
    printf,1,itime, 0,bx(i),by(i),bz(i),vx(i),vy(i),vz(i),n(i),100000.0, $
      format = '(i5, 6i3, 8f10.2)'
endfor

close,1

end

