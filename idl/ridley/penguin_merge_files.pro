
spawn, 'pwd', pwd

if (strpos(pwd,'sys_3') ge 0) then system = 'SYS3'
if (strpos(pwd,'sys_4') ge 0) then system = 'SYS4'
if (strpos(pwd,'sys_5') ge 0) then system = 'SYS5'
if (strpos(pwd,'sys_6') ge 0) then system = 'SYS6'

;Sys_3:
;   2010-12-25:
;   Latitude:       -89.984777 deg
;   Longitude:      -154.726070 deg
;
;Sys_4:
;   2011-01-21:
;   Latitude:       -89.997795 deg
;   Longitude:      123.505243 deg

systems = ['SYS3','SYS4']
lats = [-89.984777, -89.997795]
lons = ([-154.726070, 123.505243] + 360.0) mod 360.0

filelist = findfile('*.csv.gz')
nFiles = n_elements(filelist)

spawn,'wc Done',wc
wc = fix(wc)
wc = wc(0)
if (wc gt 0) then begin
   done = strarr(wc)
   close,1
   openr,1,'Done'
   line = ''
   for i=0,wc-1 do begin
      readf,1,line
      done(i) = strmid(line,0,4)+'_'+strmid(line,4,2)+'_'+strmid(line,6,2)
   endfor
   close,1
endif else done = ['nada']

IsFirst = 1
for iFile = 0,nFiles-1 do begin

   filename = filelist(iFile)
   IsFound = 0
   for i=0,wc-1 do begin
      if (strpos(filename,done(i)) ge 0) then IsFound = 1
   endfor

   if (not IsFound) then begin

      print, 'Reading file : ',filename
      spawn, 'gzip -d '+filename
      file = strmid(filename, 0, strlen(filename)-3)

      penguin_read_csv_new, file, Vars, time1, data1, error1

      if (IsFirst) then begin
         time = time1
         data = data1
         error = error1
         IsFirst = 0
      endif else begin
         time = [time,time1]
         data = [data,data1]
         error = [error,error1]
      endelse

      spawn, 'gzip '+file

   endif

endfor

if (n_elements(time) gt 0) then begin

   c_r_to_a, iTimeStart, min(time)
   c_r_to_a, iTimeEnd,   max(time)

   iTimeStart = [iTimeStart(0),iTimeStart(1),iTimeStart(2),0,0,0]
   iTimeEnd   = [iTimeEnd(0),iTimeEnd(1),iTimeEnd(2),23,59,59]

   c_a_to_r, iTimeStart, sTime
   c_a_to_r, iTimeEnd,   eTime

   nDays = fix((eTime - sTime)/(24.0*3600.0))+1

   print, 'nDays = ',nDays

   close,2
   openw,2,'Done', /append

   i = 0L

   daydata = fltarr(3,24L*3600L)
   daytime = dblarr(24L*3600L)

   for iDay = 0, nDays-1 do begin

      t = sTime + double(iDay)*24.0D*3600.0D
      c_r_to_a, iTime, t
      c_a_to_ymd, iTime, ymd
      fileout = ymd+system+'.MAG'
      print,'Writing file : ',fileout
      openw,1,fileout
      l = where(systems eq system)
      ii = l(0)
      printf,1, system, lats(ii), lons(ii), ' '+ymd+' GEOMAGNETIC nT  1Hz', $
             format = '(a,2f8.3,a)'
      nTimes = n_elements(time)
      for iT = 0L,24L*3600L-1 do begin
         c_r_to_a, itime, t
         c_a_to_ymdhms, itime, ymdhms
         xyz = 99999.0*[1.0,1.0,1.0]
         e = ' x'
         if (i lt nTimes-1) then begin
            while (time(i) lt t and i lt nTimes-2) do begin
               print, "Correcting : ",i, nTimes-2, $
                      time(i), t, t-time(i)
               i = i + 1L
            endwhile
            if (time(i) eq t) then begin
               xyz = reform(data(i,*))
               if (error(i) eq 0) then e = ' .'
               i = i + 1L
            endif
         endif
         daydata(*,iT) = xyz
         daytime(iT) = t
         printf, 1, ymdhms, xyz, e, format = '(a,3f10.3,a)'
         t = t + 1.0D
      endfor
      close,1

      l = where(abs(daydata) gt 70000.0,c)
      if (c gt 0) then daydata(l) = -1.0e32
      psfile = ymd+system+'.ps'
      setdevice,psfile,'l',5
      vars = ['Bx(nT)','By(nT)','Bz(nT)']
      overplots = [0,0,0]
      nVars = 3
      plot_data, daydata, daytime, nvars, vars, overplots
      closedevice

      if (iDay lt nDays-1) then printf,2,ymd

   endfor

   close,2

endif

end

