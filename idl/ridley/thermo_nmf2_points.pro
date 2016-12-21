
if (n_elements(filename) eq 0) then filename = 'CHAMP.location.ev1.dat'
filename = ask('filename to evaluate',filename)

if (n_elements(dir) eq 0) then dir = 'data'
dir = ask('directory for 3D files',dir)

print, 'Reading File : ',filename
close, /all
openr,1,filename
line = ''
readf,1,line
readf,1,line
readf,1,line

spawn, 'wc '+filename,wc
nLines = fix(wc)-3
nLines = nLines(0)

tmp = fltarr(7)
time = dblarr(nLines)
lat  = fltarr(nLines)
lon  = fltarr(nLines)

iLine = 0

while not eof(1) do begin

   readf,1,tmp
   itime = [tmp(0),1,tmp(1),tmp(2),tmp(3),tmp(4)]
   c_a_to_r, itime, rtime

   time(iLine) = rtime
   lat(iLine) = tmp(5)
   lon(iLine) = (tmp(6)+360.0) mod 360.0

   iLine = iLine + 1

endwhile

close,1

ymdhmsa_old = ''

for i=0,nLines-1 do begin

   RealTime = time(i)
   c_r_to_a, itime, RealTime

   itimeb = itime
   itimea = itime

   itimeb(5) = 0
   if (itimeb(4) le 29) then itimeb(4) = 0 $
   else itimeb(4) = 30

   itimea(5) = 0
   if (itimea(4) le 29) then itimea(4) = 30 $
   else begin
      itimea(3) = itimea(3)+1
      itimea(4) = 0
   endelse

   c_a_to_ymdhms, itimeb, ymdhmsb, /twoyear
   ymdhmsb = strmid(ymdhmsb,0,6)+'_'+strmid(ymdhmsb,6,6)

   c_a_to_ymdhms, itimea, ymdhmsa, /twoyear
   ymdhmsa = strmid(ymdhmsa,0,6)+'_'+strmid(ymdhmsa,6,6)

   if (ymdhmsa ne ymdhmsa_old) then begin
      print, ymdhmsb, '->',ymdhmsa

      if (ymdhmsb eq ymdhmsa_old) then begin
         fileBefore = FileAfter
         dataBefore = dataAfter
         iTImeBefore = iTimeAfter
      endif else begin
         fileBefore = dir+'/3DALL_t'+ymdhmsb+'.bin'
         read_thermosphere_file, fileBefore, nvars, nalts, nlats, nlons, $
                                 vars, dataBefore, rb, cb, bl_cnt, iTimeBefore, Version
      endelse
      fileafter = dir+'/3DALL_t'+ymdhmsa+'.bin'
      read_thermosphere_file, fileafter, nvars, nalts, nlats, nlons, $
                              vars, dataafter, rb, cb, bl_cnt, iTimeafter, Version
   endif 

   ymdhmsa_old = ymdhmsa

   c_a_to_r, iTimeBefore, rTimeBefore
   c_a_to_r, iTimeAfter, rTimeAfter

   tFac = (RealTime-rTimeBefore)/(rTimeAfter-rTimeBefore)

   GitmLons = reform(dataBefore(0,*,0,0))/!dtor
   GitmLats = reform(dataBefore(1,0,*,0))/!dtor
   GitmAlts = reform(dataBefore(2,0,0,*))/1000.0

   lo = where(GitmLons gt lon(i))
   lo = lo(0)-1
   loFac = (lon(i) - GitmLons(lo)) / (GitmLons(lo+1) - GitmLons(lo))

   la = where(GitmLats gt lat(i))
   la = la(0)-1
   laFac = (lat(i) - GitmLats(la)) / (GitmLats(la+1) - GitmLats(la))
   
   iVar = 3
   while (strpos(vars(iVar),'e-') eq -1) do iVar = iVar + 1

   ; Find specific point before time

   nmf2_mm = max(dataBefore(iVar,lo  ,la  ,*))/1.0e6
   nmf2_pm = max(dataBefore(iVar,lo+1,la  ,*))/1.0e6
   nmf2_mp = max(dataBefore(iVar,lo  ,la+1,*))/1.0e6
   nmf2_pp = max(dataBefore(iVar,lo+1,la+1,*))/1.0e6

   nmf2b = (1.0-loFac)*(1.0-laFac) * nmf2_mm + $
           (    loFac)*(1.0-laFac) * nmf2_pm + $
           (1.0-loFac)*(    laFac) * nmf2_mp + $
           (    loFac)*(    laFac) * nmf2_pp

   l = where(dataBefore(iVar,lo  ,la  ,*) eq max(dataBefore(iVar,lo  ,la  ,*)))
   hmf2_mm = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo+1,la  ,*) eq max(dataBefore(iVar,lo+1,la  ,*)))
   hmf2_pm = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo  ,la+1,*) eq max(dataBefore(iVar,lo  ,la+1,*)))
   hmf2_mp = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo+1,la+1,*) eq max(dataBefore(iVar,lo+1,la+1,*)))
   hmf2_pp = GitmAlts(l(0))

   hmf2b = (1.0-loFac)*(1.0-laFac) * hmf2_mm + $
           (    loFac)*(1.0-laFac) * hmf2_pm + $
           (1.0-loFac)*(    laFac) * hmf2_mp + $
           (    loFac)*(    laFac) * hmf2_pp

   l = where(dataBefore(iVar,lo  ,la  ,*) gt 0.9*max(dataBefore(iVar,lo  ,la  ,*)),c)
   hmf2b_mm = GitmAlts(l(0))
   hmf2a_mm = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo+1,la  ,*) gt 0.9*max(dataBefore(iVar,lo+1,la  ,*)),c)
   hmf2b_pm = GitmAlts(l(0))
   hmf2a_pm = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo  ,la+1,*) gt 0.9*max(dataBefore(iVar,lo  ,la+1,*)),c)
   hmf2b_mp = GitmAlts(l(0))
   hmf2a_mp = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo+1,la+1,*) gt 0.9*max(dataBefore(iVar,lo+1,la+1,*)),c)
   hmf2b_pp = GitmAlts(l(0))
   hmf2a_pp = GitmAlts(l(c-1))

   hmf2bb = (1.0-loFac)*(1.0-laFac) * hmf2b_mm + $
            (    loFac)*(1.0-laFac) * hmf2b_pm + $
            (1.0-loFac)*(    laFac) * hmf2b_mp + $
            (    loFac)*(    laFac) * hmf2b_pp

   hmf2ab = (1.0-loFac)*(1.0-laFac) * hmf2a_mm + $
           (    loFac)*(1.0-laFac) * hmf2a_pm + $
           (1.0-loFac)*(    laFac) * hmf2a_mp + $
           (    loFac)*(    laFac) * hmf2a_pp

   ; Find Specific point after time

   nmf2_mm = max(dataBefore(iVar,lo  ,la  ,*))/1.0e6
   nmf2_pm = max(dataBefore(iVar,lo+1,la  ,*))/1.0e6
   nmf2_mp = max(dataBefore(iVar,lo  ,la+1,*))/1.0e6
   nmf2_pp = max(dataBefore(iVar,lo+1,la+1,*))/1.0e6

   nmf2a = (1.0-loFac)*(1.0-laFac) * nmf2_mm + $
           (    loFac)*(1.0-laFac) * nmf2_pm + $
           (1.0-loFac)*(    laFac) * nmf2_mp + $
           (    loFac)*(    laFac) * nmf2_pp

   l = where(dataBefore(iVar,lo  ,la  ,*) eq max(dataBefore(iVar,lo  ,la  ,*)))
   hmf2_mm = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo+1,la  ,*) eq max(dataBefore(iVar,lo+1,la  ,*)))
   hmf2_pm = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo  ,la+1,*) eq max(dataBefore(iVar,lo  ,la+1,*)))
   hmf2_mp = GitmAlts(l(0))
   l = where(dataBefore(iVar,lo+1,la+1,*) eq max(dataBefore(iVar,lo+1,la+1,*)))
   hmf2_pp = GitmAlts(l(0))

   hmf2a = (1.0-loFac)*(1.0-laFac) * hmf2_mm + $
           (    loFac)*(1.0-laFac) * hmf2_pm + $
           (1.0-loFac)*(    laFac) * hmf2_mp + $
           (    loFac)*(    laFac) * hmf2_pp

   l = where(dataBefore(iVar,lo  ,la  ,*) gt 0.9*max(dataBefore(iVar,lo  ,la  ,*)),c)
   hmf2b_mm = GitmAlts(l(0))
   hmf2a_mm = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo+1,la  ,*) gt 0.9*max(dataBefore(iVar,lo+1,la  ,*)),c)
   hmf2b_pm = GitmAlts(l(0))
   hmf2a_pm = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo  ,la+1,*) gt 0.9*max(dataBefore(iVar,lo  ,la+1,*)),c)
   hmf2b_mp = GitmAlts(l(0))
   hmf2a_mp = GitmAlts(l(c-1))
   l = where(dataBefore(iVar,lo+1,la+1,*) gt 0.9*max(dataBefore(iVar,lo+1,la+1,*)),c)
   hmf2b_pp = GitmAlts(l(0))
   hmf2a_pp = GitmAlts(l(c-1))

   hmf2ba = (1.0-loFac)*(1.0-laFac) * hmf2b_mm + $
            (    loFac)*(1.0-laFac) * hmf2b_pm + $
            (1.0-loFac)*(    laFac) * hmf2b_mp + $
            (    loFac)*(    laFac) * hmf2b_pp

   hmf2aa = (1.0-loFac)*(1.0-laFac) * hmf2a_mm + $
            (    loFac)*(1.0-laFac) * hmf2a_pm + $
            (1.0-loFac)*(    laFac) * hmf2a_mp + $
            (    loFac)*(    laFac) * hmf2a_pp

   nmf2  = (1.0-tfac) * nmf2b  + tfac * nmf2a
   hmf2  = (1.0-tfac) * hmf2b  + tfac * hmf2a
   hmf2a = (1.0-tfac) * hmf2ab + tfac * hmf2aa
   hmf2b = (1.0-tfac) * hmf2bb + tfac * hmf2ba

   c_r_to_a, itime, RealTime
   jd = jday(itime(0),itime(1),itime(2))

   if (i eq 0) then begin
      c_a_to_ymd, itime, ymd
      l = strpos(filename,'.')
      fileout_n = strmid(filename,0,l)+'_'+ymd+'_nmf2.ccmc'
      fileout_h = strmid(filename,0,l)+'_'+ymd+'_hmf2.ccmc'
      print, 'Opening files for output : ',fileout_n, ' and ', fileout_h
      openw,2,fileout_n
      openw,3,fileout_h

;      printf,2, '#year doy hh mm ss lat lon nmf2'
;      printf,2, '#         hr mi s  deg deg cm-3'

;      printf,3, '#year doy hh mm ss lat lon hmf2  hmf2a  hmf2b'
;      printf,3, '#         hr mi s  deg deg km    km     km'

   endif

   printf,3, itime(0), jd, itime(3), itime(4), itime(5), $
          lat(i), lon(i), hmf2, hmf2a, hmf2b, $
          format = '(5i6,2f9.3,3f10.3)'
   printf,2, itime(0), jd, itime(3), itime(4), itime(5), $
          lat(i), lon(i), nmf2, $
          format = '(5i6,2f9.3,e12.5)'

endfor

close,2,3

end
