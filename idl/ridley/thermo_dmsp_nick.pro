
minlat = 50.0

filelist = findfile('3DALL*.bin')
file = filelist(0)

gitm_read_header, file, StartTime, nVars, Vars, $
                  nLons, nLats, nAlts, version

file = filelist(n_elements(filelist)-1)

gitm_read_header, file, EndTime, nVars, Vars, $
                  nLons, nLats, nAlts, version


gitm_read_bin_1var, file, gitmalts, gitmtime, nVars, Vars, ver,$
                    VarsToGet = [2]

; Grab the  altitude a two cells below the ghostcells
nAlts = n_elements(gitmalts(0,0,0,*))
GitmMaxAlt = gitmalts(0,0,0,nAlts-4)/1000.0

get_gitm_times, filelist, AllGitmTimes

DmspRootDir = '/raid3/Data/DMSP/'
read_dmsp_data, DmspRootDir, DmspData, starttime, endtime

nFiles = n_elements(DmspData)
psfile = 'compare_dmsp.ps'
setdevice, psfile, 'p', 5

SatList = [DmspData(0).satellite]

for iFile = 1, nFiles-1 do begin
   l = where(SatList eq DmspData(iFile).satellite, c)
   if (c eq 0) then SatList = [SatList, DmspData(iFile).satellite]
endfor

display, SatList
if (n_elements(iSat) eq 0) then iSat = 0
iSat = fix(ask('DMSP satellite to compare with',tostr(iSat)))

loc = where(DmspData(*).satellite eq SatList(iSat), count)

if (count gt 0) then begin

   i = loc(0)
   nPts = DmspData(i).nPts
   DmspTimes = DmspData(i).time(0:nPts-1)
   iDensity = DmspData(i).nDensity(0:nPts-1)
   oDensity = reform(DmspData(i).nDensity(0:nPts-1) * DmspData(i).Frac(0,0:nPts-1))
   iTemp = reform(DmspData(i).temperature(0,0:nPts-1))
   eTemp = reform(DmspData(i).temperature(1,0:nPts-1))
   quality = reform(DmspData(i).quality(0,0:nPts-1))
   lons = reform(DmspData(i).geopos(0,0:nPts-1))
   lats = reform(DmspData(i).geopos(1,0:nPts-1))
   alts = reform(DmspData(i).geopos(2,0:nPts-1))

   mlat = reform(DmspData(i).magpos(1,0:nPts-1))
   mlts = reform(DmspData(i).magpos(0,0:nPts-1))

   for iPt = 1, count-1 do begin

      i = loc(iPt)

      nPts = DmspData(i).nPts
      DmspTimes = [DmspTimes,DmspData(i).time(0:nPts-1)]
      iDensity = [iDensity,DmspData(i).nDensity(0:nPts-1)]
      oDensity = [oDensity,reform(DmspData(i).nDensity(0:nPts-1) * DmspData(i).Frac(0,0:nPts-1))]
      iTemp = [iTemp,reform(DmspData(i).temperature(0,0:nPts-1))]
      eTemp = [eTemp,reform(DmspData(i).temperature(1,0:nPts-1))]
      quality = [quality,reform(DmspData(i).quality(0,0:nPts-1))]

      lons = [lons,reform(DmspData(i).geopos(0,0:nPts-1))]
      lats = [lats,reform(DmspData(i).geopos(1,0:nPts-1))]
      alts = [alts,reform(DmspData(i).geopos(2,0:nPts-1))]

      mlat = [mlat,reform(DmspData(i).magpos(1,0:nPts-1))]
      mlts = [mlts,reform(DmspData(i).magpos(0,0:nPts-1))]

   endfor

   VarsToGet = [24]
   if (max(alts) gt GitmMaxAlt) then $
      alts = alts-max(alts) + GitmMaxAlt
   get_gitm_points, DmspTimes, lons, lats, alts, VarsToGet, GitmData

   ; Let's figure out when we are in the north polar region:

   loc = where(mlat gt 40.0, count)

   dt = DmspTimes(loc(1:count-1)) - DmspTimes(loc(0:count-2))

   loc2 = where(dt gt 3600.0, count2)

   maxrange = 90.0 - minlat
   print, count2

   if (count2 gt 2) then begin

      edges = [0, loc2, count-1]

      nPasses = n_elements(edges)-1

      for iPass = 0,nPasses-1 do begin

         iStart = edges(iPass)
         iEnd   = edges(iPass+1)

         if (iEnd-iStart gt 100.0) then begin

            l = loc(iStart+1:iEnd-1)

            range = 90.0 - mlat(l)
            theta = (mlts(l)-6.0)*!pi/12.0

            dtimes = DmspTimes(l)
            TimeAtMax = dtimes(where(range eq min(range)))
            TimeAtMax = TimeAtMax(0)
            delta = abs(AllGitmTimes-TimeAtMax)
            whichfile = where(delta eq min(delta))
            whichfile = filelist(whichfile(0))

            gitm_read_bin_1var, whichfile, onedata, onetime, nVars, Vars, ver,$
                                VarsToGet = [0,1,2,24]

            onealts = onedata(2,0,0,*)/1000.0
            iAlt = where(onealts ge GitmMaxAlt)
            iAlt = iAlt(0)-1
            onelats = onedata(1,0,*,iAlt)/!dtor
            onelons = onedata(0,*,0,iAlt)/!dtor

            oDmsp = oDensity(l)*1.0e6
            oGitm = gitmdata(l)

            d = [oDensity(loc)*1.0e6]
;            d = [oDmsp,oGitm]
            positives = where(d gt 0)
            
            determine_min_max, d(positives), mini, maxi

            plotdumb

            ppp = 3
            space = 0.05
            pos_space, ppp, space, sizes, nx = 3

            get_position, ppp, space, sizes, 0, pos
            
            plot, maxrange*[-1.0,1.0],maxrange*[-1.0,1.0], $
                  ystyle = 5, xstyle = 5, pos = pos, /noerase, /nodata

            v = reform(onedata(3,*,*,iAlt))

            c_r_to_a, itime, onetime
            ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
            utrot = ut * 15.0

            contour_circle, v, onelons+utrot, onelats, $
                    mini = mini, maxi = maxi, $
                    nLevels = 31, $
                    no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                    pos = pos, $
                    maxrange = maxrange, title = title, $
                    colorbar = colorbar ; , $
;                    latExtra = latExtra, lonExtra=lonExtra, $
;                    nocolor = nocolor, $
;                    nomaxmin = nomaxmin, thick = thick, average = average, $
;                    indicatemaxmin = indicatemaxmin



            plotmlt, maxrange, /no06

            x = range * cos(theta)
            y = range * sin(theta)

            oplot, x, y, thick = 4

            get_position, ppp, space, sizes, 1, pos1
            get_position, ppp, space, sizes, 2, pos2

            pos = pos1
            pos(0) = pos(0)+0.05
            pos(2) = pos2(2)

            stime = min(dtimes)
            etime = max(dtimes)
            time_axis, stime, etime, btr, etr, $
                       xtickname, xtitle, xtickv, xminor, xtickn

            plot, dtimes-stime, oGitm/1.0e11, yrange = [mini,maxi*1.1]/1.0e11, $
                  pos = pos, /noerase, $
                  xtickname = xtickname,			$
                  xtitle = xtitle,			$
                  xtickv = xtickv,			$
                  xminor = xminor,			$
                  xticks = xtickn,   $
                  ytitle = 'O+ Density (x10!U11!N)'


            oplot, oDmsp/1.0e11, linestyle = 2, min_val=0.0

            ;stop

         endif

      endfor

   endif

endif

closedevice
end
