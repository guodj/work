
if (n_elements(file) eq 0) then begin
  filelist = findfile('imf*.dat')
  file = filelist(0)
endif
file = ask('filename to change',file)

imf_read, file, time, mag, vel, den, temp, nPts, notes

print, notes

smooth = '15'
smooth = float(ask('number of minutes for smoothing',smooth))

p=strpos(file,'.dat')
fileout = strmid(file,0,p)+'_smooth_'+tostr(smooth)+'.dat'

mag_smooth = mag
vel_smooth = vel
den_smooth = den
temp_smooth = temp
  
smooth = smooth*60.0

for i=0,npts-1 do begin
  
   l = where(time ge time(i)-smooth/2 and time lt time(i)+smooth/2,c)

   if (c gt 0) then begin
      iStart = l(0)
      iEnd   = l(c-1)
   endif

;   if i-smooth/2 ge 0 then istart = i-smooth/2 else istart = 0
;   if i+smooth/2 le npts-1 then iend = i+smooth/2 else iend = npts-1
  
   for j=0,2 do begin
      mag_smooth(j,i) = mean(mag(j,istart:iend))
      vel_smooth(j,i) = mean(vel(j,istart:iend))
   endfor
   den_smooth(i) = mean(den(istart:iend))
   temp_smooth(i) = mean(temp(istart:iend))
  
endfor

openw,1,fileout

for i=0,n_elements(notes)-1 do printf,1,notes(i)

printf,1,''
printf,1,'Data smoothed using program imf_smooth.pro'
printf,1,'Smoothing time : ',smooth
printf,1,''
printf,1,'#START'

for i=0,npts-1 do begin

   c_r_to_a, itime, time(i)

   printf,1, format = '(i5,5i3,i4,3f8.2,4f9.2,f13.1)', $
          itime, 0, mag_smooth(*,i), $
          vel_smooth(*,i), den_smooth(i), temp_smooth(i)

endfor

close,1

end
