
filelist = findfile('3DALL*.bin')

nfiles = n_elements(filelist)

ff = 528
lf = nfiles-1
;lf = ff + 100
;lf = 300

iVar = [0,1,2,16,17]

iE_ = 3
iN_ = iE_+1

nErrors = 180
AngleErrors = findgen(nErrors)

NorthErrors = fltarr(nErrors)
EastErrors = fltarr(nErrors)

meanmag = 0.0

for ifile = ff, lf do begin

   file = filelist(iFile)

   print, 'Reading file : ',file

   gitm_read_bin_1var,file, data, time, nVars, Vars, version, VarsToGet = iVar

   if (iFile eq ff) then begin
      nlons = n_elements(data(0,*,0,0))
      nlats = n_elements(data(1,0,*,0))
      nalts = n_elements(data(2,0,0,*))
      lats = reform(data(1,*,*,0))/!dtor
      lons = reform(data(0,*,*,0))/!dtor
      Alts = reform(data(2,0,0,*))/1000.0
      l = where(Alts gt 245.0)
      iAlt = l(0)
      Alt = Alts(iAlt)
      iVar = [16,17]

      loc = where(lats gt -80.0 and lats lt -60.0)

   endif

   East = data(iE_,*,*,iAlt)
   East = East(loc)
   North = data(iN_,*,*,iAlt)
   North = North(loc)
   mag = sqrt(East^2 + North^2)
   meanmag = meanmag + mean(mag)
   angle = acos(East/mag)
   aLoc = where(North lt 0, c)
   if (c gt 0) then angle(aLoc) = 2*!pi - angle(aLoc)

   for iA = 0,nErrors-1 do begin
      e = randomn(s,nLons*nLats)
      NewAngle = angle + e*AngleErrors(iA)*!dtor
      NewEast = mag * cos(NewAngle)
      NewNorth = mag * sin(NewAngle)
      NorthError = sqrt(mean((North - NewNorth)^2))
      EastError = sqrt(mean((East - NewEast)^2))
      NorthErrors(iA) = NorthErrors(iA) + NorthError
      EastErrors(iA) = EastErrors(iA) + EastError
   endfor

   iE_ = 0
   iN_ = iE_+1

endfor

NorthErrors = NorthErrors / (lf-ff+1)
EastErrors = EastErrors / (lf-ff+1)
mag = sqrt(NorthErrors^2+EastErrors^2)
meanmag = meanmag/(lf-ff+1)

setdevice, 'wind_errors.ps','p',5

plot, AngleErrors, mag, xtitle = 'Pointing Error (Degrees)', $
      ytitle = 'Wind Error (m/s)', thick = 4, xrange = [0,20], $
      pos = [0.1, 0.3, 0.99, 0.8]
xyouts,  AngleErrors(15), mag(15)+1, 'Mean Total Error', alignment = 1.0

oplot, AngleErrors, NorthErrors, linestyle = 2, thick = 2
xyouts,  AngleErrors(15), NorthErrors(15)+1, 'Mean Northward Error', charsize = 0.8, orient = 25.0
oplot, AngleErrors, EastErrors, linestyle = 3, thick = 2
xyouts,  AngleErrors(15), EastErrors(15)-2, 'Mean Eastward Error', charsize = 0.8, orient = 21.0

oplot, AngleErrors, AngleErrors*0.0+meanmag/10.0, linestyle = 1
xyouts, 19, meanmag/10.0+0.5, 'Desired Maximum Error', alignment = 1.0

oplot, [3.0,3.0], [0,100], linestyle = 1
xyouts, 2.9, 40.0, 'Attitude Pointing Capability', orient = 90.0

closedevice

end
