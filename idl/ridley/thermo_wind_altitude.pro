filelist = findfile('3DALL*.bin')

nfiles = n_elements(filelist)

ff = 82
lf = nfiles-1
lf = 82

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
;      lats = reform(data(1,1:nLons-2,1:nLats-2,0)) ; /!dtor
;      lons = reform(data(0,1:nLons-2,1:nLats-2,0)) ; /!dtor
      lats = reform(data(1,*,*,0)) ; /!dtor
      lons = reform(data(0,*,*,0)) ; /!dtor
      Alts = reform(data(2,0,0,*))/1000.0
      l = where(Alts gt 245.0)
      iAlt1 = l(0)
      Alt1 = Alts(iAlt1)

      l = where(Alts gt 400.0)
      iAlt2 = l(0)
      Alt2 = Alts(iAlt2)

      iVar = [16,17]
   endif

;   East1 = reform(data(iE_,1:nLons-2,1:nLats-2,iAlt1))
;   East2 = reform(data(iE_,1:nLons-2,1:nLats-2,iAlt2))
   East1 = reform(data(iE_,*,*,iAlt1))
   East2 = reform(data(iE_,*,*,iAlt2))

;   North1 = reform(data(iN_,1:nLons-2,1:nLats-2,iAlt1))
;   North2 = reform(data(iN_,1:nLons-2,1:nLats-2,iAlt2))
   North1 = reform(data(iN_,*,*,iAlt1))
   North2 = reform(data(iN_,*,*,iAlt2))

   mag1 = sqrt(East1^2 + North1^2)
   mag2 = sqrt(East2^2 + North2^2)

   psfile = 'test.ps'

;   value = (mag2-mag1)/mag2*100.0
   value = mag2-mag1

   mini = -100.0
   maxi = 100.0
   title = 'Neutral Wind ('+tostr(Alt2)+'km-'+tostr(Alt1)+'km)'
   colortitle = 'Magnitude Difference (m/s)'
   vn = North2-North1
   ve = East2-East1

   maxrange = 40.0

   thermo_threeplot, psfile, value, time, lons, lats, mini, maxi, $
                      title, colortitle, maxrange, $
                      vn = vn, ve = ve

   iE_ = 0
   iN_ = iE_+1

endfor

end

