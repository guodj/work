
function get_length, file

  close,1
  openr,1,file

  line = ''
  nVars = -1
  while (strpos(line,'elements') lt 0 and not eof(1)) do readf,1,line
  if (strpos(line,'elements') lt 0) then iLength = -1 $
  else begin
     l = strpos(line,'=')+1
     iLength = long(strmid(line,l,10))

     readf,1,line  ; dead line
     readf,1,line  ; header line

     vars = strsplit(line,' ',/extract)
     nVars = n_elements(vars)

  endelse

  result = [iLength,nVars]
  return, result

end

pro read_log_file, file, time, data, vars

  result = get_length(file)

  nLines = result(0)
  nVars = result(1)

  data = fltarr(nVars,nLines)
  time = dblarr(nLines)

  close,1
  openr,1,file

  line = ''
  while (strpos(line,'elements') lt 0 and not eof(1)) do readf,1,line
  if (strpos(line,'elements') gt 0) then begin
     readf,1,line  ; dead line
     readf,1,line  ; header line
     vars = strsplit(line,' ',/extract)
     readf,1,data
     for iT=0L,nLines-1 do begin
        iTime = fix(data(1:6,iT))
        c_a_to_r, iTime, rTime
        time(iT) = rTime
     endfor
  endif


end

; look at ../../_db.php

dbfile = '../../_db.php'
spawn, 'pwd',pwd
pwd = pwd[0]
p = strpos(pwd,'/sats/')
satellite = strmid(pwd,p+6,strlen(pwd))
p = strpos(satellite,'/')
satellite = strmid(satellite,0,p)

factor = 1.0
if strpos(satellite,'Geotail') ge 0 then factor = 0.1

spawn, 'grep '+satellite+' '+dbfile,results
CdfVar = results(3)
p = strpos(CdfVar,'=>')
CdfVar = strmid(CdfVar,p+4,strlen(CdfVar))
p = strpos(CdfVar,"'")
CdfVar = strmid(CdfVar,0,p)

filelist_cdf = findfile('*.cdf')
nFilesCdf = n_elements(filelist_cdf)

nPtsCdf = 0

for iFile = 0,nFilesCdf-1 do begin

   cdffile = filelist_cdf(iFile)
   id = cdf_open(cdffile)
   info = cdf_inquire(id)

   natt = info.natts
   nvar = info.nvars
   if nvar eq 0 then begin
      nvar = info.nzvars
      zvar = 1
   endif else zvar = 0

   nrows = info.maxrec

   if (nrows lt 1) then begin
      cdf_control, id, variable=0, /zvariable,get_var_info=v
      nrows = v.maxrec
      if (nrows le 0) then nrows = v.maxrecs
   endif

   nPtsCdf = nPtsCdf + nrows

   cdf_close, id

endfor

AllCdfData = fltarr(3,nPtsCdf)
AllCdfTimes = dblarr(nPtsCdf)

n = 0

for iFile = 0,nFilesCdf-1 do begin

   cdffile = filelist_cdf(iFile)
   print, 'reading cdf file : ', cdffile
   id = cdf_open(cdffile)
   info = cdf_inquire(id)

   natt = info.natts
   nvar = info.nvars
   if nvar eq 0 then begin
      nvar = info.nzvars
      zvar = 1
   endif else zvar = 0

   nrows = info.maxrec

   if (nrows lt 1) then begin
      cdf_control, id, variable=0, /zvariable,get_var_info=v
      nrows = v.maxrec
      if (nrows le 0) then nrows = v.maxrecs
   endif

   vars = strarr(nvar)

   epoch  = -1
   epoch0 = -1
   for i=0,nvar-1 do begin

      if (zvar eq 0) then r = cdf_varinq(id,i)            $
      else r = cdf_varinq(id,i,/zvariable)
      vars(i) = r.name
      tmp = mklower(vars(i))
      ; stupid themis calls its time variable time, but it is offset
      ; from the variable epoch0, so we need to read in both...
      if ((strpos(tmp,'epoch') gt -1 or strpos(tmp,'time') gt -1) and epoch lt 0) then epoch = i
      if (strpos(tmp,'epoch0') gt -1 and epoch0 lt 0) then epoch0 = i

   endfor

   cdf_varget, id, vars(epoch), times, rec_count=nrows

   swap = 0
   if (epoch0 gt -1) then begin
      cdf_varget, id, vars(epoch0), epoch0time
      if (max(epoch0time) lt 1.0e10) then swap = 1
   endif
   if (max(times) lt 1.0e10) then swap = 1

   cdf_varget, id, CdfVar, CdfData, rec_count=nrows

   if (swap eq 1) then begin
      CdfTime = swap_endian(times)
      CdfData = swap_endian(CdfData)
      if (epoch0 gt -1) then begin
         CdfTime = swap_endian(epoch0time) + CdfTime*1000.0
         ;print, cdftime
      endif
   endif else CdfTime = times

   cdf_close,id

   s = size(CdfTime)
   if (s(0) eq 2) then CdfTime = reform(CdfTime(0,*))

   for iRow = 0, nRows-1 do begin
      cdf_epoch, CdfTime(iRow), year, month, day, hour, minute, second, milli, /break
      iTime = [ year, month, day, hour, minute, second]
      ;print, itime
      ;stop
      c_a_to_r, iTime, rtime
      CdfTime(iRow) = rtime
   endfor

   AllCdfData(*,n:n+nrows-1) = CdfData*factor
   AllCdfTimes(n:n+nrows-1) = CdfTime

   n = n + nrows

endfor

if (strpos(mklower(satellite),'cluster') eq -1) then begin
   filelist_txt = findfile('*GSM*.txt')
endif else begin
   filelist_txt = findfile('*GSE*.txt')
   print, 'CLUSTER! Looking for GSE file!!!'
endelse

nFilesTxt = n_elements(filelist_txt)

iLenMax = -1
for iFile = 0,nFilesTxt-1 do begin
   result = get_length(filelist_txt(iFile))
   iLen = result(0)
   nVars = result(1)
   if (iLen gt iLenMax) then iLenMax = iLen
   print, 'reading file ',filelist_txt(iFile),' with values ',iLen,nVars
endfor

AllData = fltarr(nFilesTxt, nVars, iLenMax)
AllTime = dblarr(nFilesTxt, iLenMax)

for iFile = 0,nFilesTxt-1 do begin
   read_log_file, filelist_txt(iFile), time, data, vars
   nLines = n_elements(data(0,*))
   AllTime(iFile,0:nLines-1) = time
   AllData(iFile,*,0:nLines-1) = data
endfor

l = where(alltime gt 0)

c_r_to_a, itime_start, min(alltime(l))
c_r_to_a, itime_end, max(alltime(l))

c_a_to_s, itime_start, stime_start
c_a_to_s, itime_end, stime_end

print, 'Start time : ',stime_start
print, 'end time : ',stime_end

if (n_elements(start_time) eq 0) then start_time = strmid(stime_start,0,15)
start_time = ask('start time of plotting',start_time)

if (strlen(start_time) lt 9) then $
  start_time = strmid(stime_start,0,9)+' '+start_time

;--------------------------------------------------------------
; I got sick of typing in the ending date, so if the date is
; the same, assume the user just wants to enter the time
;--------------------------------------------------------------

sdate = strmid(start_time,0,9)
if strpos(stime_end,sdate) gt -1 then 					$
  end_time_default = strmid(stime_end,10,5)				$
else end_time_default = strmid(stime_end,0,15)

if (n_elements(end_time) eq 0) then end_time = end_time_default
end_time   = ask('end time of plotting',end_time)

;--------------------------------------------------------------
; If the user entered a short string, assume it is just a time
; and add the date on the front
;--------------------------------------------------------------

if (strlen(end_time) lt 9) then end_time = strmid(start_time,0,9)+' '+end_time

;--------------------------------------------------------------
; Now figure out where in the file these things are, with the
; default to give the user everything
;--------------------------------------------------------------

c_s_to_a, itime_start, start_time
c_a_to_r, itime_start, stime

c_s_to_a, itime_end, end_time
c_a_to_r, itime_end, etime

if (etime lt stime) then begin
   mid = (stime+etime)/2.0
   stime = mid - 12.0*3600.0
   etime = stime + 24.0*3600.0
endif

display, vars
if (n_elements(iVar) eq 0) then iVar = 15
iVar = ask('which variable to compare',tostr(iVar))

possible_vars = ['bx','by','bz']

for i=0,2 do begin
   if (strpos(mklower(vars(iVar)),possible_vars(i)) ge 0) then iVarCdf = i
endfor

; -------------------
; compute rms differences
; -------------------

nPtsCdf = n_elements(AllCdfTimes)
delta  = fltarr(nFilesTxt, nPtsCdf)-1.0e32
deltaN = fltarr(nFilesTxt, nPtsCdf)-1.0e32

rms  = fltarr(nFilesTxt)
ave  = fltarr(nFilesTxt)
Nrms = fltarr(nFilesTxt)
nptsRms = lonarr(nFilesTxt)

overlapdata  = fltarr(nFilesTxt, nPtsCdf)-1.0e32
overlapmodel = fltarr(nFilesTxt, nPtsCdf)-1.0e32

for iTime = 0L, nPtsCdf-1 do begin

   if (AllCdfData(iVarCDF,iTime) gt -1.0e31) then begin

      for iFile = 0, nFilesTxt-1 do begin

         d = abs(AllCdfTimes(iTime)-alltime(iFile,*))
         l = where(d eq min(d))
         if (d(l(0)) lt 60.0 and alldata(iFile,iVar,l(0)) gt -1.0e31) then begin
            delta(iFile,iTime) = alldata(iFile,iVar,l(0)) - AllCdfData(iVarCdf,iTime)
            overlapdata(iFile,iTime) = AllCdfData(iVarCdf,iTime)
            overlapmodel(iFile,iTime) = alldata(iFile,iVar,l(0))
            if (AllCdfData(iVarCdf,iTime) ne 0.0) then begin
               deltaN(iFile,iTime) = $
                  100.0*(alldata(iFile,iVar,l(0)) - AllCdfData(iVarCdf,iTime))/$
                  AllCdfData(iVarCdf,iTime)
            endif

         endif

      endfor

   endif

endfor

v = reform(alldata(*,iVar,*))

l = where(overlapdata gt -1.0e31,c)
if (c gt 0) then begin
   combined = [overlapdata(l), overlapmodel(l)]
   if (n_elements(mini) eq 0) then determine_min_max, combined, mini, maxi
endif else begin
   mini = -1.0
   maxi = 1.0
endelse
mini = float(ask('min value to plot',string(mini)))
maxi = float(ask('max value to plot',string(maxi)))

for iFile = 0, nFilesTxt-1 do begin
   rms(iFile) = -1.0
   ave(iFile) = -1.0
   Nrms(iFile) = -1.0
   d = reform(delta(iFile,*))
   l = where(d gt -1.0e31,c)
   nPtsRms(iFile) = c
   if (c gt 0) then begin
      rms(iFile) = sqrt(median(d(l)^2))
      ave(iFile) = median(d(l))
   endif
   d = reform(deltaN(iFile,*))
   l = where(d gt -1.0e31,c)
   if (c gt 0) then Nrms(iFile) = sqrt(median(d(l)^2))
endfor

c_r_to_a, itime, stime
c_a_to_ymdhms, itime, sTimeString

c_r_to_a, itime, etime
c_a_to_ymdhms, itime, eTimeString

l = min([strlen(vars(iVar)),2])
psfile = satellite+'_'+strmid(vars(iVar),0,l)+'_'+sTimeString+'_'+eTimeString+'.ps'
difffile = satellite+'_'+strmid(vars(iVar),0,l)+'_'+sTimeString+'_'+eTimeString+'.txt'
setdevice,psfile,'p',5

makect,'bry'

dColor = 245.0/(nFilesTxt+1)

close,2
openw,2,difffile

for iFile = 0,nFilesTxt-1 do begin

   if (iFile eq 0) then begin

      time_axis, stime, etime, btr, etr, xtickname, xtitle, $
                 xtickv, xminor, xtickn

      l = where(alltime(iFile,*) ge stime and alltime(iFile,*) le etime, c)

      s = size(v)
      if (s(0) eq 2) then vv = reform(v(iFile,*)) else vv = v

      plot, alltime(iFile,l)-stime, vv(l), $
            xrange = [btr,etr], xstyle = 1, $
            yrange = [mini,maxi], $
            ytitle = vars(iVar), $
            xtickname = xtickname,			$
            xtitle = xtitle,			$
            xtickv = xtickv,			$
            xminor = xminor,			$
            xticks = xtickn, /nodata, pos = [0.1,0.45,0.95,0.95], $
            title = satellite+' '+vars(iVar)

   endif
   l = where(alltime(iFile,*) ge stime and alltime(iFile,*) le etime, c)
   if (c gt 0) then begin
      c = dColor + float(iFile)*dColor
      oplot, alltime(iFile,l)-stime, vv(l), color = c, thick = 2
      dx = (etr-btr)/10.0
      dy = (maxi-mini)/nFilesTxt*0.8
      plots, [btr+0.1*dx, btr+dx], [mini,mini]-dy*(iFile+2), color = c, thick = 4
      s = filelist_txt(iFile)+' '+string(rms(iFile),format='(f8.2)')+'nT'+$
          ', '+string(ave(iFile),format='(f9.2)')+'nT'+ $
          ' ('+string(Nrms(iFile),format='(f8.2)')+'%) - '+ $
          tostr(nPtsRms(iFile),5)+' points'
      xyouts, btr+dx*1.1, mini-dy*(iFile+2), s
      printf,2,s
   endif
endfor

close,2

; plot CDF data

oplot, AllCdfTimes-stime, AllCdfData(iVarCDF,*), min_val = mini, $
       thick = 5, color = 0

closedevice

end
