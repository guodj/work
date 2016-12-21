filename = 'imf_by-.dat'

bx_mean = 0.0 
bx_amp  = 0.0

by_mean = -5.0
by_amp  =  0.0

bz_mean =    0.0
bz_amp  =    0.0

vx_mean = -400.0
vx_amp  =    0.0

vy_mean =   0.0
vy_amp  =   0.0

vz_mean =   0.0
vz_amp  =   0.0

n_mean =  5.0
n_amp  =  0.0

itime = [2011, 12, 21, 0, 0, 0]
c_a_to_r, itime, basetime

dt = 15.0 * 60.0

nDays  = 4
nTimes = nDays * (24.0*3600.0)/dt + 1
nFrequencies = 15

MaxFrequency = 60.0
MinFrequency =  4.0

frequency = fltarr(nTimes)
phase     = fltarr(nTimes)

time = findgen(nTimes) * dt / 60.0

nT = nTimes/nFrequencies
for i = 0, nFrequencies-2 do begin
    frequency(i*nT:(i+1)*nT-1) = $
      i*(MaxFrequency-MinFrequency)/(nFrequencies-1) + MinFrequency
endfor
i = nFrequencies-1
frequency(i*nT:nTimes-1) = $
  i*(MaxFrequency-MinFrequency)/(nFrequencies-1) + MinFrequency

wave = (time / frequency) * 2*!pi

for i = 1, nFrequencies-2 do begin
    print, 'before : ',i*nT,(i+1)*nT-1, wave(i*nt-2:i*nT+2)
    wave(i*nT:(i+1)*nT-1) = wave(i*nT:(i+1)*nT-1) - $
      wave(i*nT) + wave(i*nT-1)
    print, 'after : ',i*nT,(i+1)*nT-1, wave(i*nt-2:i*nT+2)
endfor

i = nFrequencies-1
print, 'before : ',i*nT,(i+1)*nT-1, wave(i*nt-2:i*nT+2)
wave(i*nT:nTimes-1) = wave(i*nT:nTimes-1) - $
  wave(i*nT) + wave(i*nT-1)
print, 'after : ',i*nT,(i+1)*nT-1, wave(i*nt-2:i*nT+2)

bx = bx_mean + bx_amp * sin(wave)
by = by_mean + by_amp * sin(wave)
bz = bz_mean + bz_amp * sin(wave)
vx = vx_mean + vx_amp * sin(wave)
vy = vy_mean + vy_amp * sin(wave)
vz = vz_mean + vz_amp * sin(wave)
n  =  n_mean +  n_amp * sin(wave)

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

