
if n_elements(filein) eq 0 then begin
   filelist = findfile('power*')
   filein = filelist(0)
endif

filein = ask('file to read an rearrange',filein)

if n_elements(fileout) eq 0 then fileout = 'power_normalized.dat'
fileout = ask('output file',fileout)

if n_elements(dt) eq 0 then dt = 5.0
dt = float(ask('time (dt) to normalize (minutes)',string(dt)))

if n_elements(Combine) eq 0 then Combine = 0
Combine = fix(ask('whether you want to combine North & South (yes=1)',$
                  tostr(Combine)))

power_read, filein, time, power, North, nPts, notes

stime = min(time)
c_r_to_a, itime, stime
itime([3,4,5]) = 0
c_a_to_r, itime, starttime

etime = max(time)
c_r_to_a, itime, etime
itime([3,4,5]) = [23,59,59]
c_a_to_r, itime, endtime

nTimes = fix((endtime-starttime)/(dt*60.0))

timeN = dblarr(nTimes)
timeS = dblarr(nTimes)

PowerN = fltarr(nTimes)
PowerS = fltarr(nTimes)

l = where(North eq 1,nN)
time0N = time(l)
power0N = power(l)

l = where(North eq 0,nS)
time0S = time(l)
power0S = power(l)

time  = dblarr(nTimes*2)
power = fltarr(nTimes*2)
North = intarr(nTimes*2)

for i=0,nTimes-1 do begin

   timeN(i) = starttime + double(i)*dt*60.0
   timeS(i) = timeN(i) + dt*60.0/2.0

   d = timen(i) - time0n

   print, timen(i), min(abs(d))

   if (timeN(i) lt min(time0N)) then PowerN(i) = Power0N(0) $
   else if (timeN(i) gt max(time0N)) then PowerN(i) = Power0N(nN-1) $
   else begin
      l = where(time0N ge timeN(i)-dt/2 and time0N lt timeN(i)+dt/2, c)
      if (c gt 0) then PowerN(i) = mean(PowerN(l)) $
      else begin
         l = where(time0N ge timeN(i))
         iL = l(0)-1
         iU = l(0)
         x = (timeN(i) - time0N(iL)) / (time0N(iU) - time0N(iL))
         PowerN(i) = (1.0-x)*Power0N(iL) + x*PowerN(iU)
      endelse
   endelse

   if (timeS(i) lt min(time0S)) then PowerS(i) = Power0S(0) $
   else if (timeS(i) gt max(time0S)) then PowerS(i) = Power0S(nS-1) $
   else begin
      l = where(time0S ge timeS(i)-dt/2 and time0S lt timeS(i)+dt/2, c)
      if (c gt 0) then PowerS(i) = mean(PowerS(l)) $
      else begin
         l = where(time0S ge timeS(i))
         iL = l(0)-1
         iU = l(0)
         x = (timeS(i) - time0S(iL)) / (time0S(iU) - time0S(iL))
         PowerS(i) = (1.0-x)*Power0S(iL) + x*PowerS(iU)
      endelse
   endelse

   time(i*2) = timeN(i)
   power(i*2) = PowerN(i)
   North(i*2) = 1

   time(i*2+1) = timeS(i)
   power(i*2+1) = PowerS(i)
   North(i*2+1) = 0

   print, Combine
   if (Combine) then begin
      power(i*2) = (PowerN(i) + PowerS(i))/2.0
      power(i*2+1) = (PowerN(i) + PowerS(i))/2.0
   endif

endfor

power_write, fileout, time, power, North, notes

end
