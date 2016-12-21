
file = 'imf.dat'

imf_read, file, time, mag, vel, den, temp, nPts, notes

nTimes = n_elements(time)

dt = time(1:nTimes-1)-time(0:nTimes-2)
meandt = mean(dt)

OneHour = 3600.0

nHours = fix((max(time)-min(time))/OneHour+0.5)

sTime = min(time) + OneHour/2.0
eTime = max(time) - OneHour/2.0
nHours = fix((eTime-sTime)/OneHour)

print, 'Mean dt : ',meandt

for i=0,2 do print, 'std b : ',i,stdev(mag(i,*))
for i=0,2 do print, 'std v : ',i,stdev(vel(i,*))
print, 'Mean den : ',stdev(den)
print, 'Mean temp : ',stdev(temp)

magstd = fltarr(3,nHours)
magmean = fltarr(3,nHours)
velstd = fltarr(3,nHours)
denstd = fltarr(nHours)
tempstd = fltarr(nHours)

for iHour= 0,nHours-1 do begin

   t = sTime + iHour*OneHour
   l = where(time gt t-OneHour/2 and time le t+OneHour/2, c)
   if (c gt 1) then begin
      for i = 0,2 do begin
         magstd(i,iHour) = stdev(mag(i,l))
         magmean(i,iHour) = mean(mag(i,l))
         velstd(i,iHour) = stdev(vel(i,l))
      endfor
      tempstd(iHour) = stdev(temp(l))
      denstd(iHour) = stdev(den(l))
   endif

endfor

for i=0,2 do print, 'hourly std b : ',i,mean(magstd(i,*))
for i=0,2 do print, 'hourly std v : ',i,mean(velstd(i,*))
print, 'hourly std n : ',mean(denstd(where(denstd gt 0)))
print, 'hourly std t : ',mean(tempstd(where(tempstd gt 0)))

end



