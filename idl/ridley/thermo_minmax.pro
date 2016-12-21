
filelist = findfile('3DALL*bin')
nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

    filename = filelist(iFile)

    print, 'Reading file ',filename

    read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
      vars, data, rb, cb, bl_cnt, iTime, Version

    if (iFile eq 0) then begin

       minmax = fltarr(3,nVars,nFiles,nAlts-4)
       lats = reform(data(1,*,*,0))
       lons = reform(data(0,*,*,0))
       nLons = n_elements(lons(*,0))
       nLats = n_elements(lats(0,*))
       dlon = lons
       dlon(0:nlons-2,*) = lons(1:nLons-1,*) - lons(0:nLons-2,*)
       dlon(nlons-1,*) = dlon(0,*)

       dlat = lats
       dlat(*,1:nlats-2) = (lats(*,2:nLats-1) - lats(*,0:nLats-3))/2.0
       dlat(*,0) = lats(*,1) - lats(*,0)
       dlat(*,nLats-1) = lats(*,nLats-1) - lats(*,nLats-2)
       area = 6372000.0*6372000.0 * cos(lats) * dlon * dlat

    endif

    for iVar = 0, nVars-1 do begin

       for iAlt = 2, nAlts-3 do begin

          minmax(0,iVar,iFile,iAlt-2) = min(data(iVar,*,*,iAlt))
          minmax(1,iVar,iFile,iAlt-2) = max(data(iVar,*,*,iAlt))
          minmax(2,iVar,iFile,iAlt-2) = mean(data(iVar,*,*,iAlt)*area)/mean(area)

       endfor

    endfor

endfor

openw,1,'.minmax'
printf,1,nVars,nAlts-4

for iVar = 0, nVars-1 do begin

   printf,1,vars(iVar)
   for iAlt = 2, nAlts-3 do begin
      mi = min( minmax(0,iVar,*,iAlt-2))
      ma = max( minmax(1,iVar,*,iAlt-2))
      me = mean(minmax(2,iVar,*,iAlt-2))
      printf,1,mi,ma,me
   endfor

endfor

close,1

end



