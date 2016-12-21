
pro line_plot, direction, input, maxwind, GitmLatsAll, GitmUTs, $
               gitmwind, sdiwind, gitme

  nFiles = n_elements(GitmUTs)

  setdevice, direction+'_'+input+'_sdi_line.ps','p',5

  space = 0.1
  ppp = 1
  pos_space, ppp, space, sizes

  mini = min(GitmLatsAll)
  maxi = max(GitmLatsAll)

  l=where(GitmLatsAll gt mini)
  dl = GitmLatsAll(l(0)) - mini

  mini = mini-dl
  maxi = maxi+dl

  r = maxwind/dl

  get_position, ppp, space, sizes, 0, pos

  stime = min(GitmUTs)
  etime = max(GitmUTs)

  time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

  makect,'wb'

  ; Figure out how many lats there are:

  nLats = 0
   la = GitmLatsAll
  while (min(la) lt 1000.0) do begin
     l = where(la eq min(la))
     la(l) = 1000.0
     nLats++
  endwhile

  eden = fltarr(nFiles, nLats)
  x = fltarr(nFiles, nLats)
  y = fltarr(nFiles, nLats)

  for i=0,nLats-1 do x(*,i) = GitmUTs-stime

  iLat = 0
  la = GitmLatsAll
  while (min(la) lt 1000.0) do begin
     l = where(la eq min(la))
     y(*,iLat) = min(la)
     for iFile=0,nFiles-1 do eden(iFile,iLat) = mean(Gitme(iFile,l))
     la(l) = 1000.0
     iLat++
  endwhile

  mine = 0.0
  maxe = max(gitme)
  range = maxe-mine
  maxe = maxe + r/2
  range = maxe-mine
  levels = findgen(31)/30 * range + mine

  l = where(eden lt levels(2),c)
  if (c gt 0) then eden(l) = levels(2)

  contour, eden, x, y, xstyle = 5, ystyle = 5, $
           pos = pos, levels = levels, /fill

  ctpos = pos
  ctpos(1) = ctpos(1) - 0.07
  ctpos(3) = ctpos(1) + 0.02

  title = 'Electron Density (/m!U3!N)'
  ncolors = 255
  minmax = mm(levels)
  plotct, ncolors, ctpos, minmax, title, /bottom

  plot, [btr, etr], [mini, maxi], $
        ystyle = 5, xstyle = 1, $
        xtickname = xtickname,			$
        xtitle = xtitle,			$
        xtickv = xtickv,			$
        xminor = xminor,			$
        xticks = xtickn,   $
        pos = pos, /nodata, $
        title = direction+' -> '+input, /noerase

  xyouts, -etr/20.0, (mini+maxi)/2.0, 'Latitude (deg)', $
          orient = 90, align = 0.5, charsize = 1.1

  la = GitmLatsAll

  iLat = 0

  makect,'bry'

  while (min(la) lt 1000.0) do begin

     l = where(la eq min(la))

     gitmave = fltarr(nFiles)
     sdiave = fltarr(nFiles)

     for iFile=0,nFiles-1 do begin

        gitmave(iFile) = mean(GitmWind(iFile,l))
        sdiave(iFile) = mean(SdiWind(iFile,l))

     endfor

     gitmave = gitmave/r
     sdiave  = sdiave/r

     offset = mini + iLat*dl + dl

     oplot, gitmUTs-stime, gitmave+offset, color = 10, thick = 4
     oplot, gitmUTs-stime, sdiave+offset,  color = 200, thick = 4

     oplot, [btr,etr], [offset,offset], linestyle = 1, thick=2

     sLat = string(min(la),format='(f5.2)')
     xyouts, -etr/100.0, offset, sLat, orient = 90, align = 0.5

     la(l) = 1000.0

     iLat++

  endwhile

  dx = etr/100.0

  offset1 = mini + iLat*dl + dl
  iLat--
  offset2 = mini + iLat*dl + dl

  plots, [etr+dx,etr+dx], [offset1, offset2], thick = 2
  xyouts, etr+dx, offset2-dl/10.0, $
          tostr(maxwind)+' (m/s)', orient = 270

  r = (maxi-mini)/10.0
  offset1 = (maxi+mini)/2.0
  offset2 = offset1 - r
  plots, [etr+dx,etr+dx], [offset1, offset2], thick = 2, color = 10
  xyouts, etr+dx*2, offset2, 'GITM'

  r = (maxi-mini)/10.0
  offset1 = (maxi+mini)/2.0 - r*3
  offset2 = offset1 - r
  plots, [etr+dx,etr+dx], [offset1, offset2], thick = 2, color = 200
  xyouts, etr+dx*2, offset2, 'SDI'

  closedevice

end

input = 'weimer'
version = '01'
Dir = 'run.'+input+version

restore, 'Ridley_Winds_20121124.sav'

lons = lons+360.0

MinLon = min(lons)
MaxLon = max(lons)

MinLat = min(lats)
MaxLat = max(lats)

filelist = findfile(Dir+'/*.bin')
nFiles = n_elements(filelist)

VarsToGet = [0,1,2,16,17,18,33]

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, 'Reading ',file
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet

   ut = (time mod (24.0*3600.0))/3600.0

   if (iFile eq 0) then begin

      GitmLons = reform(data(0,*,*,0))/!dtor
      GitmLats = reform(data(1,*,*,0))/!dtor
      GitmAlts = reform(data(2,0,0,*))/1000.0

      GitmUTs  = dblarr(nFiles)

      list = where(GitmLons gt MinLon and $
                   GitmLons lt MaxLon and $
                   GitmLats gt MinLat and $
                   GitmLats lt MaxLat, nPts)

      print, nPts

      SdiLoc = intarr(nPts)

      SdiVeAll = fltarr(nFiles, nPts)
      SdiVnAll = fltarr(nFiles, nPts)
      SdiVzAll = fltarr(nFiles, nPts)

      SdiLatsAll = fltarr(nPts)
      SdiLonsAll = fltarr(nPts)

      GitmVeAll = fltarr(nFiles, nPts)
      GitmVnAll = fltarr(nFiles, nPts)
      GitmVzAll = fltarr(nFiles, nPts)
      GitmeAll = fltarr(nFiles, nPts)

      GitmLatsAll = fltarr(nPts)
      GitmLonsAll = fltarr(nPts)

      for iPt = 0, nPts-1 do begin

         dlo = (GitmLons(list(iPt)) - lons)*cos(lats*!dtor)
         dla = GitmLats(list(iPt)) - lats

         d = dla^2 + dlo^2

         l = where(d eq min(d))
         SdiLoc(iPt) = l(0)

         GitmLatsAll(iPt) = GitmLats(list(iPt))
         GitmLonsAll(iPt) = GitmLons(list(iPt))

         SdiLatsAll(iPt) = lats(SdiLoc(iPt))
         SdiLonsAll(iPt) = lons(SdiLoc(iPt))

      endfor

      dAlt = abs(240.0-GitmAlts)
      l = where(dAlt eq min(dAlt))
      iAlt = l(0)

   endif

   GitmVe = reform(data(3,*,*,iAlt))
   GitmVn = reform(data(4,*,*,iAlt))
   GitmVz = reform(data(5,*,*,iAlt))
   Gitme  = reform(data(6,*,*,iAlt))

   GitmUTs(iFile) = time

   dt = abs(ut - UT_TIME_DECIMAL_HOURS)
   l = where(dt eq min(dt))
   iT = l(0)

   print, ut, UT_TIME_DECIMAL_HOURS(iT)

   for iPt = 0, nPts-1 do begin

      j = SdiLoc(iPt)

      GitmVeAll(iFile, iPt) = GitmVe(list(iPt))
      GitmVnAll(iFile, iPt) = GitmVn(list(iPt))
      GitmVzAll(iFile, iPt) = GitmVz(list(iPt))
      GitmeAll(iFile, iPt)  = Gitme(list(iPt))
      
      SdiVeAll(iFile, iPt) = winds(SdiLoc(iPt),iT).zonal
      SdiVnAll(iFile, iPt) = winds(SdiLoc(iPt),iT).meridional
      SdiVzAll(iFile, iPt) = winds(SdiLoc(iPt),iT).vertical

   endfor

endfor

maxwind = 500.0

line_plot, 'North', input, maxwind, GitmLatsAll, GitmUTs, $
               GitmVnAll, SdiVnAll, GitmeAll

line_plot, 'East', input, maxwind, GitmLatsAll, GitmUTs, $
               GitmVeAll, SdiVeAll, GitmeAll

line_plot, 'Vertical', input, 50.0, GitmLatsAll, GitmUTs, $
               GitmVzAll, SdiVzAll, GitmeAll

stop


end
