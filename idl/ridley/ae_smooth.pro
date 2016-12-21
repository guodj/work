
if (n_elements(file) eq 0) then begin
  filelist = findfile('sme*.dat')
  file = filelist(0)
endif
file = ask('filename to change',file)

ae_read, file, time, ae

c_r_to_a, itime_start, min(time)
c_a_to_s, itime_start, stime_start

c_r_to_a, itime_end, max(time)
c_a_to_s, itime_end, stime_end

if (n_elements(start_date) eq 0) then start_date = strmid(stime_start,0,9)
if (n_elements(end_date) eq 0) then end_date = strmid(stime_end,0,9)

print,"start time in file : ",stime_start
print,"end time in file   : ",stime_end

start_date = ask('starting date',start_date)
end_date   = ask('ending date',end_date)

c_s_to_a, itime, start_date
c_a_to_r, itime, starttime

c_s_to_a, itime, end_date+' 23:59:59'
c_a_to_r, itime, endtime

if (n_elements(starttime) gt 0) then begin
   l = where(time ge starttime,c)
   if (c gt 0) then begin
      time = time(l)
      ae = ae(l,*)
   endif
endif

if (n_elements(endtime) gt 0) then begin
   l = where(time le endtime,c)
   if (c gt 0) then begin
      time = time(l)
      ae = ae(l,*)
   endif
endif

if (n_elements(smooth) eq 0) then smooth = 15
smooth = float(ask('number of minutes for smoothing',tostr(smooth)))

smooths = smooth*60.0

ae_temp = ae

nPts = n_elements(time)

for i=0L,npts-1 do begin
  
   l = where(time ge time(i)-smooths/2 and time lt time(i)+smooths/2,c)

   if (c gt 0) then begin
      iStart = l(0)
      iEnd   = l(c-1)
   endif

   for j=0,7 do ae(i,j) = mean(ae_temp(l,j))

endfor

p=strpos(file,'.dat')
fileout = strmid(file,0,p)+'_smooth_'+tostr(smooth)+'.dat'

close,2
openw,2,fileout
for i=0L,nPts-1 do begin
   c_r_to_a, itime, time(i)
   printf,2,itime,ae(i,*), format = '(6i5,3f8.0,f8.1,f10.0,3f10.1)'
endfor
close,2

end
