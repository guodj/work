
pro get_fpi_vel, FpiAz, FpiEl, FpiTime, FpiVelocity, FpiTemp, $
                 FlagV, FlagT, GitmE, GitmN, GitmT, GitmTime, $
                 El0, dEl, Az0, dAz, OutVel, OutTemp, OutTime, $
                 OutFlagV, OutFlagT, $
                 OutGitmV, OutGitmT, OutGitmTime

  l = where(abs(FpiEl) lt 5.0, nZenith)
  if (nZenith gt 0) then begin
     zenith = FpiVelocity(l)
     zenithtime = FpiTime(l)
  endif

  FpiAzAlt = FpiAz
  if (Az0 lt dAz) then begin
     l = where(FpiAz gt 360+Az0-dAz,c)
     if (c gt 0) then fpiAzAlt(l) = fpiAzAlt(l)+360
  endif

  l = where(fpiAzAlt gt Az0-dAz and fpiAzAlt lt Az0+dAz and $
            fpiEl gt El0-dEl and fpiEl lt El0+dEl ,c)

  if (c gt 0) then begin

     los = fpiVelocity(l)
     PointTime = fpitime(l)
     ze = fpiEl(l)*!dtor
     if (nZenith gt 0) then $
        w = interpolate_mine(PointTime, zenith, zenithTime) $
     else w = 0.0
     si = 1.0
     if (Az0 gt 135 and Az0 lt 315) then si = -1.0
     OutVel = si * (los - w*cos(ze))/sin(ze)
     OutTemp = FpiTemp(l)
     OutTime = FpiTime(l)
     OutFlagV = FlagV(l)
     OutFlagT = FlagT(l)

     OutGitmT = GitmT(l)
     OutGitmTime = GitmTime(l)
     az = FpiAz(l)*!dtor
     OutGitmV = si * (GitmN(l) * cos(az) + GitmE(l) * sin(az))

  endif else begin

     OutVel = [0.0]
     OutTemp = [0.0]
     OutTime = [0.0]
     OutFlagV = [0]
     OutFlagT = [0]
     OutGitmT = [0.0]
     OutGitmV = [0.0]
     OutGitmTime = [0.0]

  endelse

end


pro get_var_at_loc, varIn, iLon, iLat, xLon, xLat, varOut

  varOut = (1.0-xLon) * (1.0-xLat) * varIn(iLon,  iLat  ) + $
           (    xLon) * (1.0-xLat) * varIn(iLon+1,iLat  ) + $
           (1.0-xLon) * (    xLat) * varIn(iLon,  iLat+1) + $
           (    xLon) * (    xLat) * varIn(iLon+1,iLat+1)
  
end



filelist = findfile('minime*.txt');
nFiles = n_elements(filelist)

iStats = intarr(nFiles)

i = 0

iStatSave = [fix(strmid(filelist(i),6,2))]
for i=1,nFiles-1 do begin
   iStatTry = fix(strmid(filelist(i),6,2))
   l = where(iStatSave eq iStatTry,c)
   if (c eq 0) then iStatSave = [iStatSave,iStatTry]
endfor

nStats = n_elements(iStatSave)

FPINums = [1,2,3,4,5,6,7,8,9]
FPINames = ['car','caj','unk','unk','uao','par','eku','ann','vti']

FPILats =       [-7.28, -6.88, 0.0, 0.0, 40.1097, 35.2417, 37.7447, 42.40, 37.23]
FPILons = 360.0-[36.52, 38.56, 0.0, 0.0, 88.2042, 82.7286, 84.2936, 83.91, 80.42]

if (n_elements(gitmdir) eq 0) then gitmdir = ''
gitmdir = ask('gitm directory (data.xxx or nothing for data)',gitmdir)
if (strlen(gitmdir) gt 0) then dir = 'data.'+gitmdir else dir = 'data'

for iS = 0, nStats-1 do begin

   iStat = iStatSave(iS)
   l = where(FPINums eq iStat,c)
   if (c eq 0) then begin
      print, 'Station '+tostr(iStat,2)+' not found! Skipping!'
   endif else begin

      fpiloc = [FPILats(iStat-1), FPILons(iStat-1)]

      filelist = findfile('minime0'+tostr(FPINums(iStat-1))+'*.txt')
      display, filelist

      fpi_read_makela_all, filelist, fpi, allfpidata

      FpiAz = AllFpiData.Azimuth
      FpiEl = AllFpiData.Zenith
      FpiTime = AllFpiData.Time
      FpiVelocity = AllFpiData.Velocity
      FpiTemperature = AllFpiData.Temperature
      FpiFlagV = AllFpiData.FlagV
      FpiFlagT = AllFpiData.FlagT

      l = where(FpiEl gt 60.0, c)
      if (c gt 0) then FpiEl(l) = 0.0

      sb = 6372.0
      sa = sb + 250.0

      A = (180.0 - FpiEl)*!dtor
      B = asin(sin(A)*sb/sa)
      C = !pi - (A+B)

      LatList = fpiloc(0)+ (C/!dtor) * cos(FpiAz*!dtor)
      LonList = fpiloc(1)+ (C/!dtor) * sin(FpiAz*!dtor)/cos(LatList*!dtor)

      CardEl = [45.0, 45.0, 45.0, 45.0, 0.0]
      CardAz = [0.0, 90.0, 180.0, 270.0, 0.0]

      A = (180.0 - CardEl)*!dtor
      B = asin(sin(A)*sb/sa)
      C = !pi - (A+B)

      LatCard = fpiloc(0)+ (C/!dtor) * cos(CardAz*!dtor)
      LonCard = fpiloc(1)+ (C/!dtor) * sin(CardAz*!dtor)/cos(LatCard*!dtor)

      filelist = findfile(dir+'/3DALL*bin')

      nPts = n_elements(FpiTime)

      gitm_read_header, filelist, gitmtimes, nVars, Vars, $
                        nLons, nLats, nAlts, version

      loc = where(gitmtimes ge min(FpiTime) and $
                  gitmtimes lt max(FpiTime)+3600.0, count)

      if count gt 0 then begin

         iVar = [0,1,2,15,16,17,18]
         gitmUn   = fltarr(nPts)
         gitmUe   = fltarr(nPts)
         gitmUu   = fltarr(nPts)
         gitmTn   = fltarr(nPts)
         gitmTime = dblarr(nPts)

         gitmUnCard = fltarr(5,count)
         gitmUeCard = fltarr(5,count)
         gitmTnCard = fltarr(5,count)

         iPt = 0

         gitmtimes = dblarr(count)

         for iFile = loc(0), loc(count-1) do begin

            file = filelist(iFile)
            print, file
            gitm_read_bin_1var,file, data, time1, nVars, $
                               Vars, version, VarsToGet = iVar

            if (iFile eq loc(0)) then begin

               lons = reform(data(0,*,0,0)/!dtor)
               dlon = lons(1)-lons(0)
               iLon = (LonList-lons(0))/dlon
               xLon = iLon - fix(iLon)
               iLon = fix(iLon)

               iLonC = (LonCard-lons(0))/dlon
               xLonC = iLonC - fix(iLonC)
               iLonC = fix(iLonC)

               lats = reform(data(1,0,*,0)/!dtor)
               dlat = lats(1)-lats(0)
               iLat = (LatList-lats(0))/dlat
               xLat = iLat - fix(iLat)
               iLat = fix(iLat)

               iLatC = (LatCard-lats(0))/dlat
               xLatC = iLatC - fix(iLatC)
               iLatC = fix(iLatC)

               alts = reform(data(2,0,0,*)/1000.0)
               l = where(alts gt 240.0)
               iAlt = l(0)

            endif

            Tn1 = reform(data(3,*,*,iAlt)) 
            Ue1 = reform(data(4,*,*,iAlt))
            Un1 = reform(data(5,*,*,iAlt))
            Uu1 = reform(data(6,*,*,iAlt))

            gitmtimes(iFile-loc(0)) = time1
            for i=0,4 do begin
               get_var_at_loc, Un1, iLonC(i), iLatC(i), $
                               xLonC(i), xLatC(i), Out
               gitmUnCard(i,iFile-loc(0)) = Out
               get_var_at_loc, Ue1, iLonC(i), iLatC(i), $
                               xLonC(i), xLatC(i), Out
               gitmUeCard(i,iFile-loc(0)) = Out
               get_var_at_loc, Tn1, iLonC(i), iLatC(i), $
                               xLonC(i), xLatC(i), Out
               gitmTnCard(i,iFile-loc(0)) = Out
            endfor

            if (iFile gt loc(0)) then begin

               l = where(FpiTime ge time0 and FpiTime lt time1, c)

               if (c gt 0) then begin

                  for ii = 0, c-1 do begin

                     li = l(ii)
                     t = FpiTime(li)
                     xt = (t-time0)/(time1-time0)

                     get_var_at_loc, Un0, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UnOut0
                     get_var_at_loc, Un1, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UnOut1
                     UnInt = (1-xt)*UnOut0 + xt*UnOut1

                     get_var_at_loc, Ue0, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UeOut0
                     get_var_at_loc, Ue1, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UeOut1
                     UeInt = (1-xt)*UeOut0 + xt*UeOut1

                     get_var_at_loc, Uu0, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UuOut0
                     get_var_at_loc, Uu1, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), UuOut1
                     UuInt = (1-xt)*UuOut0 + xt*UuOut1

                     get_var_at_loc, Tn0, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), TnOut0
                     get_var_at_loc, Tn1, iLon(li), iLat(li), $
                                     xLon(li), xLat(li), TnOut1
                     TnInt = (1-xt)*TnOut0 + xt*TnOut1

                     gitmUn(li) = UnInt
                     gitmUe(li) = UeInt
                     gitmUu(li)   = UuInt
                     gitmTn(li)   = TnInt
                     gitmTime(li) = t
                     iPt++

                  endfor

               endif

            endif

            Tn0 = Tn1
            Ue0 = Ue1
            Un0 = Un1
            Uu0 = Uu1
            time0 = time1

         endfor

      endif

      pos3 = [0.1, 0.05, 1.0, 0.325]
      pos2 = [0.1, 0.35, 1.0, 0.625]
      pos1 = [0.1, 0.65, 1.0, 0.925]

      mt = min(fpitime)
      c_r_to_a, itime, mt
      itime(3:5) = 0.0
      c_a_to_r, itime, stime
      c_a_to_s, itime, StrDate
      StrDate = strmid(StrDate,0,9)


      mt = max(fpitime)
      c_r_to_a, itime, mt
      itime(3) = 23
      itime(4) = 59
      itime(5) = 59
      c_a_to_r, itime, etime

      stime = min(gitmtimes)
      etime = max(gitmtimes)

      nDayStart = double(long(stime/3600.0/24.0))*3600.0*24.0
      nDayEnd  = double(long(etime/3600.0/24.0)+1)*3600.0*24.0

;      t = (fpitime - stime)

      FpiAz = AllFpiData.Azimuth
      FpiEl = AllFpiData.Zenith
      FpiTime = AllFpiData.Time
      FpiVelocity = AllFpiData.Velocity
      FpiTemperature = AllFpiData.Temperature
      FpiFlagV = AllFpiData.FlagV
      FpiFlagT = AllFpiData.FlagT

      get_fpi_vel, FpiAz, FpiEl, FpiTime, FpiVelocity, FpiTemperature, $
                   FpiFlagV, FpiFlagT, GitmUe, GitmUn, GitmTn, GitmTime, $
                   45.0, 15.0, 0.0, 30.0, $
                   FpiNorth, FpiNorthTemp, FpiNorthTime, $
                   FpiNorthVFlag, FpiNorthTFlag, $
                   GitmNorth, GitmNorthTemp, GitmNorthTime

      get_fpi_vel, FpiAz, FpiEl, FpiTime, FpiVelocity, FpiTemperature, $
                   FpiFlagV, FpiFlagT, GitmUe, GitmUn, GitmTn, GitmTime, $
                   45.0, 15.0, 90.0, 30.0, $
                   FpiEast, FpiEastTemp, FpiEastTime, $
                   FpiEastVFlag, FpiEastTFlag, $
                   GitmEast, GitmEastTemp, GitmEastTime

      get_fpi_vel, FpiAz, FpiEl, FpiTime, FpiVelocity, FpiTemperature, $
                   FpiFlagV, FpiFlagT, GitmUe, GitmUn, GitmTn, GitmTime, $
                   45.0, 15.0, 180.0, 30.0, $
                   FpiSouth, FpiSouthTemp, FpiSouthTime, $
                   FpiSouthVFlag, FpiSouthTFlag, $
                   GitmSouth, GitmSouthTemp, GitmSouthTime

      get_fpi_vel, FpiAz, FpiEl, FpiTime, FpiVelocity, FpiTemperature, $
                   FpiFlagV, FpiFlagT, GitmUe, GitmUn, GitmTn, GitmTime, $
                   45.0, 15.0, 270.0, 30.0, $
                   FpiWest, FpiWestTemp, FpiWestTime, $
                   FpiWestVFlag, FpiWestTFlag, $
                   GitmWest, GitmWestTemp, GitmWestTime


      fpiV = [fpinorth,fpisouth,fpieast,fpiwest]
      l = where(abs(fpiV) lt 1000.0,c)
      if (c gt 0) then fpiV = fpiV(l)

      determine_min_max, abs(fpiV), mini, maxi_f, percent = 95
      determine_min_max, abs([GitmUn,GitmUe]), mini, maxi_g, percent = 95

      ymax = max([maxi_f,maxi_g])
      ymax = float(floor(ymax/100.0)+1)*100.0
      ymin = -ymax
      dy = ymax-ymin

      meant = mean(fpitime)
      c_r_to_a, itime, meant
      c_a_to_ymd, itime, ymd

      print, FPINames(iStat-1)+'_winds_'+ymd+'.ps'
      setdevice, FPINames(iStat-1)+'_winds_'+ymd+'.ps','p',5

      makect,'mid'

;      stime = min(gitmtimes)
;      etime = max(gitmtimes)

      nt = (fpinorthtime - stime)
      st = (fpisouthtime - stime)

      time_axis, stime, etime, btr, etr, xtickname, xtitle, $
                 xtickv, xminor, xtickn

      xmin = btr
      xmax = etr
      dx = xmax-xmin

      plot, [btr,etr], [ymin,ymax], /nodata, $
            xtitle = ' ', ytitle = 'North Velocity (m/s)', $
            xstyle = 1, yrange = [ymin,ymax], ystyle = 1, $
            thick = 4, pos = pos1, $
	  xtickname = xtickname,			$
	  xtickv = xtickv,			$
	  xminor = xminor,			$
	  xticks = xtickn,   $            
            title = 'FPI Data From '+mkupper(FPINames(iStat-1))+'; GITM Run '+mkupper(gitmdir)

      l = where(fpinorthVflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, nt(l), fpinorth(l), 3600.0, $
                                        linestyle = 0, thick = 4

      l = where(fpinorthVflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, nt(l), fpinorth(l), 3600.0, $
                                        linestyle = 0, thick = 2, color = 240

      l = where(fpisouthVflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, st(l), fpisouth(l), 3600.0, $
                                        linestyle = 2, thick = 4
      l = where(fpisouthVflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, st(l), fpisouth(l), 3600.0, $
                                        linestyle = 2, thick = 2, color = 240

      oplot, mm(t), [0.0,0.0], linestyle = 1

      l = where(gitmtimes gt stime)
      tg = (gitmtimes - stime)
      oplot, tg, gitmUnCard(0,l), thick = 4, color = 10 
      oplot, tg, gitmUnCard(2,l), thick = 4, color = 10, linestyle = 2 

;      l = where(gitmnorthtime gt stime)
;      tgn = (gitmnorthtime - stime)/3600.0
;      oplot, tgn(l), gitmnorth(l), thick = 4, color = 10, psym = 4
;      l = where(gitmsouthtime gt stime)
;      tgs = (gitmsouthtime - stime)/3600.0
;      oplot, tgs(l), gitmsouth(l), thick = 4, color = 10, psym = 5

      l = where(gitmnorthtime gt stime and fpinorthVflag le 1,c)
      if (c gt 0) then $
         rmsnorth = sqrt(mean((gitmnorth(l)-fpinorth(l))^2)) $
      else rmsnorth = 0.0

      l = where(gitmsouthtime gt stime and fpisouthVflag le 1,c)
      if (c gt 0) then $
         rmssouth = sqrt(mean((gitmsouth(l)-fpisouth(l))^2)) $
      else rmssouth = 0.0

      xs = xmin + dx*0.05
      xe = xs + dx*0.1
      ys = ymin+dy*0.12
      oplot, [xs,xe], [ys,ys], thick = 4
      xyouts, xe+dx*0.01, ys, mkupper(FPINames(iStat-1))+' North'

      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4, linestyle = 2
      xyouts, xe+dx*0.01, ys, mkupper(FPINames(iStat-1))+' South'

      xs = xmin + dx*0.3
      xe = xs + dx*0.1
      ys = ymin+dy*0.12
      oplot, [xs,xe], [ys,ys], thick = 4, color = 10
      str = ' (RMS: '+string(rmsnorth,format='(f6.1)')+' m/s)'
      xyouts, xe+dx*0.01, ys, 'GITM North'+str, color = 10

      xs = xmin + dx*0.3
      xe = xs + dx*0.1
      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4, color = 10, linestyle = 2
      str = ' (RMS: '+string(rmssouth,format='(f6.1)')+' m/s)'
      xyouts, xe+dx*0.01, ys, 'GITM South'+str, color = 10

      for i=nDayStart-stime, nDayEnd-stime, 24*3600.0 do begin
         oplot, [i,i], [ymin, ymax], thick = 3, linestyle = 1
;         c_r_to_a, it, stime + double(i)*3600.0
;         c_a_to_s, it, sd
;         sd = strmid(sd,0,9)
;         xyouts, i+dx*0.01, ymax-dy*0.15, sd
      endfor
      

;------------------------------------------------------------

      et = (fpieasttime - stime)
      wt = (fpiwesttime - stime)

      plot, [btr,etr], [ymin,ymax], /nodata, $
            xtitle = ' ', ytitle = 'East Velocity (m/s)', $
            xstyle = 1, yrange = [ymin,ymax], ystyle = 1, $
            thick = 4, pos = pos2, /noerase, $
	  xtickname = xtickname,			$
	  xtickv = xtickv,			$
	  xminor = xminor,			$
	  xticks = xtickn

      l = where(fpieastVflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, et(l), fpieast(l), 3600.0, $
                                        linestyle = 0, thick = 4
      l = where(fpieastVflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, et(l), fpieast(l), 3600.0, $
                                        linestyle = 0, thick = 2, color = 240

      l = where(fpiwestVflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, wt(l), fpiwest(l), 3600.0, $
                                        linestyle = 2, thick = 4
      l = where(fpiwestVflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, wt(l), fpiwest(l), 3600.0, $
                                        linestyle = 2, thick = 2, color = 240

      l = where(gitmtimes gt stime)
      tg = (gitmtimes - stime)
      oplot, tg, gitmUeCard(1,l), thick = 4, color = 10 
      oplot, tg, gitmUeCard(3,l), thick = 4, color = 10, linestyle = 2 

;      l = where(gitmnorthtime gt stime)
;      tgn = (gitmnorthtime - stime)/3600.0
;      oplot, tgn(l), gitmnorth(l), thick = 4, color = 10, psym = 4
;      l = where(gitmsouthtime gt stime)
;      tgs = (gitmsouthtime - stime)/3600.0
;      oplot, tgs(l), gitmsouth(l), thick = 4, color = 10, psym = 5

      l = where(gitmeasttime gt stime and fpieastVflag le 1,c)
      if (c gt 0) then $
         rmseast = sqrt(mean((gitmeast(l)-fpieast(l))^2)) $
      else rmseast = '0.0'

      l = where(gitmwesttime gt stime and fpiwestVflag le 1,c)
      if (c gt 0) then $
         rmswest = sqrt(mean((gitmwest(l)-fpiwest(l))^2)) $
      else rmswest = '0.0'

      oplot, mm(t), [0.0,0.0], linestyle = 1

      xs = xmin + dx*0.05
      xe = xs + dx*0.1
      ys = ymin+dy*0.12
      oplot, [xs,xe], [ys,ys], thick = 4
      xyouts, xe+dx*0.01, ys, mkupper(FPINames(iStat-1))+' East'

      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4, linestyle = 2
      xyouts, xe+dx*0.01, ys, mkupper(FPINames(iStat-1))+' West'

      xs = xmin + dx*0.3
      xe = xs + dx*0.1
      ys = ymin+dy*0.12
      oplot, [xs,xe], [ys,ys], thick = 4, color = 10
      str = ' (RMS: '+string(rmseast,format='(f6.1)')+' m/s)'
      xyouts, xe+dx*0.01, ys, 'GITM East'+str, color = 10

      xs = xmin + dx*0.3
      xe = xs + dx*0.1
      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4, color = 10, linestyle = 2
      str = ' (RMS: '+string(rmswest,format='(f6.1)')+' m/s)'
      xyouts, xe+dx*0.01, ys, 'GITM West'+str, color = 10

      for i=nDayStart-stime, nDayEnd-stime, 24*3600.0 do begin
         oplot, [i,i], [ymin, ymax], thick = 3, linestyle = 1
;         c_r_to_a, it, stime + double(i)*3600.0
;         c_a_to_s, it, sd
;         sd = strmid(sd,0,9)
;         xyouts, i+dx*0.01, ymax-dy*0.15, sd
      endfor

;------------------------------------------------------------

      fpiT = [fpinorthTemp,fpisouthTemp,fpieastTemp,fpiwestTemp]
      meanFpiT = median(fpiT)

      determine_min_max, fpiT, mini_f, maxi_f, percent = 95
      determine_min_max, gitmTnCard, mini_g, maxi_g, percent = 95

      ymin = min([mini_f,mini_g])*0.8
      ymax = max([maxi_f,maxi_g])*1.1

      plot, [btr,etr], [ymin,ymax], /nodata, $
            ytitle = 'Temperature (K)', $
            xrange = mm(t), xstyle = 1, yrange = [ymin,ymax], ystyle = 1, $
	  xtickname = xtickname,			$
	  xtickv = xtickv,			$
	  xminor = xminor,			$
	  xticks = xtickn, $
            xtitle = xtitle,$
            thick = 4, pos = pos3, /noerase

      l = where(fpinorthTflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, nt(l), fpinorthTemp(l), 3600.0, thick = 4
      l = where(fpinorthTflag eq 1, c)
      if (c gt 1) then $
         oplot_with_gaps, nt(l), fpinorthTemp(l), 3600.0, thick = 2, color = 240

      l = where(fpisouthTflag lt 2, c)
      if (c gt 1) then $
         oplot_with_gaps, st(l), fpisouthTemp(l), 3600.0, linestyle = 1, thick = 4
      l = where(fpisouthTflag eq 1, c)
      if (c gt 1) then $
         oplot_with_gaps, st(l), fpisouthTemp(l), 3600.0, linestyle = 1, thick = 2, color = 240

      l = where(fpieastTflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, et(l), fpieastTemp(l), 3600.0, linestyle = 2, thick = 4
      l = where(fpieastTflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, et(l), fpieastTemp(l), 3600.0, linestyle = 2, thick = 2, color = 240

      l = where(fpiwestTflag lt 2, c)
      if (c gt 1) then oplot_with_gaps, wt(l), fpiwestTemp(l), 3600.0, linestyle = 3, thick = 4
      l = where(fpiwestTflag eq 1, c)
      if (c gt 1) then oplot_with_gaps, wt(l), fpiwestTemp(l), 3600.0, linestyle = 3, thick = 2, color = 240

;zt = (fpi.ZenithTime - stime)/3600.0
;oplot, wt, fpi.zenithTemp, linestyle = 4, thick = 4

      oplot, tg, gitmTnCard(0,*), thick = 4, color = 10 
      oplot, tg, gitmTnCard(1,*), thick = 4, color = 10, linestyle = 1 
      oplot, tg, gitmTnCard(2,*), thick = 4, color = 10, linestyle = 2 
      oplot, tg, gitmTnCard(3,*), thick = 4, color = 10, linestyle = 3 
      oplot, tg, gitmTnCard(4,*), thick = 4, color = 10, linestyle = 4 

      xs = xmin + dx*0.05
      xe = xs + dx*0.1
      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4
      xyouts, xe+dx*0.01, ys, mkupper(FPINames(iStat-1))+' Temps'


      l = where(gitmnorthtime gt stime and fpinorthTflag lt 1,c)
      if (c gt 0) then $
         diffTnorth = (gitmnorthtemp(l)-fpinorthtemp(l))^2 $
      else diffTnorth = 0.0

      l = where(gitmeasttime gt stime and fpieastTflag lt 1,c)
      if (c gt 0) then $
         diffTeast = (gitmeasttemp(l)-fpieasttemp(l))^2 $
      else diffTeast = 0.0

      l = where(gitmsouthtime gt stime and fpisouthTflag lt 1,c)
      if (c gt 0) then $
         diffTsouth = (gitmsouthtemp(l)-fpisouthtemp(l))^2 $
      else diffTsouth = 0.0

      l = where(gitmwesttime gt stime and fpiwestTflag lt 1,c)
      if (c gt 0) then $
         diffTwest = (gitmwesttemp(l)-fpiwesttemp(l))^2 $
      else diffTwest = 0.0

      rms = sqrt(mean([diffTnorth,diffTeast,diffTsouth,diffTwest]))
      nrms = rms/meanFpiT*100.0
      str = ' (RMS: '+string(rms,format='(f7.1)')+' K / '+string(nrms,format='(f7.1)')+'%)'

      xs = xmin + dx*0.35
      xe = xs + dx*0.1
      ys = ymin+dy*0.05
      oplot, [xs,xe], [ys,ys], thick = 4, color = 10
      xyouts, xe+dx*0.01, ys, 'GITM Temps'+str, color = 10

      for i=nDayStart-stime, nDayEnd-stime, 24*3600.0 do begin
         oplot, [i,i], [ymin, ymax], thick = 3, linestyle = 1
;         c_r_to_a, it, stime + double(i)*3600.0
;         c_a_to_s, it, sd
;         sd = strmid(sd,0,9)
;         xyouts, i+dx*0.01, ymax-dy*0.15, sd
      endfor

      closedevice

   endelse

endfor

end
