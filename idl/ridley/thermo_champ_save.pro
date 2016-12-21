ChampDir = '/bigdisk1/Data6/Data/CHAMP/'

reread = 1
if (n_elements(GitmData) gt 0) then begin
    answer = ask('whether to re-read gitm data','n')
    if (strpos(mklower(answer),'n') gt -1) then reread = 0
endif

if (reread) then begin

   filelist = file_search('cham_*.bin')
   gitm_read_bin, filelist, GitmData, GitmTime, nVars, Vars, version, /skip
   nTimes = n_elements(GitmTime)
   GitmAlts = reform(GitmData(0,2,0,0,*))/1000.0
   nAlts = n_elements(GitmAlts)
   GitmRho = reform(GitmData(*,3,0,0,*))

endif

GitmLons = reform(GitmData(*,0,0,0,0))
GitmLats = reform(GitmData(*,1,0,0,0))

c_r_to_a, itime, GitmTime(0)
c_r_to_a, ietime, max(GitmTime)
itime(3:5) = 0
ietime(3:5) = [23,59,59]

c_a_to_r,  itime, stime
c_a_to_r, ietime, etime

nDays = round((etime-stime)/3600.0/24.0)
doy0 = jday(itime(0), itime(1), itime(2))
iYear = itime(0)
sYear = tostr(iYear)

t = GitmTime-stime
hour = t/3600.0
GitmLocalTime = (GitmLons/!dtor/15.0 + hour) mod 24.0

; read champ data

for cday = 0, ndays - 1 do begin

   doy = doy0 + cday
   champfile = champdir+sYear+'/Density_3deg_'+$
               tostr((iYear mod 100),2)+'_'+ $
               chopr('00'+tostr(doy),3)+'.ascii'

   print, 'Reading '+champfile
   read_champ_density, champfile, ChampTime, ChampData, ChampVars, note

   if (n_elements(champtime) gt 0) then begin
   
      if (cDay eq 0) then begin
         ChampAlt = reform(champdata( 6,*))
         ChampRho = reform(champdata(11,*))
         ChampErr = reform(champdata(15,*))
         ChampTim = ChampTime
      endif else begin
         ChampAlt = [ChampAlt,reform(champdata( 6,*))]
         ChampRho = [ChampRho,reform(champdata(11,*))]
         ChampErr = [ChampErr,reform(champdata(15,*))]
         ChampTim = [ChampTim,ChampTime]
      endelse

   endif

endfor

ChampTime = ChampTim

; interpolar on to gitm times

ChampIntAlt = fltarr(nTimes)
ChampIntRho = fltarr(nTimes)
ChampIntErr = fltarr(nTimes)
GitmIntRho  = fltarr(nTimes)

for itime = 0L, nTimes-1 do begin

   dt = GitmTime(iTime) - ChampTime
   loc = where(dt lt 0, c)
   champdt = 0.0
   if (c eq 0) then begin
           ; haven't reached the start of the file yet
      i = 1
      x = 1.0
   endif else begin
      if (GitmTime(iTime) gt max(ChampTime)) then begin
         i = n_elements(ChampTime)
         x = 0.0
      endif else begin
         i = loc(0)
         if (i eq 0) then i=1
         x = (champtime(i) - GitmTime(itime)) / (champtime(i)-champtime(i-1))
      endelse
   endelse

   ChampIntRho(iTime)   = $
      ((1-x) * ChampRho(i) + x*ChampRho(i-1))*1.0e12
   ChampIntAlt(iTime)  = $
      (1-x) * ChampAlt(i) + x * ChampAlt(i-1)
   ChampIntErr(iTime)      = $
      (1-x) * ChampErr(i) + x * ChampErr(i-1)

   loc = where(GitmAlts gt ChampIntAlt(iTime),c)
   if (c gt 0) then begin

      i = loc(0)
      x = (ChampIntAlt(iTime) - GitmAlts(i-1)) / $
          (GitmAlts(i) - GitmAlts(i-1))

      GitmIntRho(iTime) = $
         exp((1.0 - x) * alog(GitmRho(iTime,i-1)) + $
             (      x) * alog(GitmRho(iTIme,i)))*1.0e12

   endif

endfor

stime = min(GitmTime)
etime = max(GitmTime)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

c_r_to_a, itime, stime
c_a_to_ymd, itime, ymd

save, file = 'gitm_champ'+ymd+'.save', $
      GitmIntRho, GitmTime, ChampIntRho, ChampIntAlt, ChampIntErr, $
      GitmLons, GitmLats, GitmLocalTime

end
