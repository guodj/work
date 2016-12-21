
pro despike_mags, data, time, lat, quality

  missing = -1.0e32

  lats = findgen(91)
  maxb = 500.0 + 1500.0*exp(-(lats-68)^2/100.0)
  maxdbdt = 2.0 + 8.0*exp(-(lats-70)^2/100.0)

  nPts  = n_elements(data(0,*,0))
  nDays = n_elements(data(0,0,*))
  nPtsTotal = nPts*nDays

  b    = fltarr(3, nPtsTotal)
  t    = dblarr(   nPtsTotal)
  dbdt = fltarr(3, nPtsTotal)
  for iComp = 0,2 do begin
      for iDay = 0, nDays-1 do begin
          b(iComp, iDay*nPts:(iDay+1)*nPts-1) = data(iComp,*,iDay)
          if (iComp eq 0) then $
            t(iDay*nPts:(iDay+1)*nPts-1) = time(*,iDay)
      endfor
      l = where(b(iComp,*) gt -1.0e31,count)
      if (count gt 0) then $
        dbdt(iComp,l(0:count-2)) = $
        (b(iComp, l(1:count-1)) - b(iComp, l(0:count-2))) / $
        (t(       l(1:count-1)) - t(       l(0:count-2)))
  endfor

  quality = intarr(nPtsTotal)

  for iComp = 0,2 do begin
      loc_good = where(b(iComp,*) gt -1.0e31,count)
      if (count gt 0) then begin

          b(iComp,loc_good) = b(iComp,loc_good) - median(b(iComp,loc_good))

          loc_bad = where(abs(b(iComp,loc_good)) gt maxb(abs(lat)), count_bad)
          if (count_bad gt 0) then begin
              print, count_bad," Bad Point(s) found in Maximum check"
              quality(loc_good(loc_bad)) = 2
          endif

          loc_bad = where(abs(dbdt(iComp,loc_good)) gt maxdbdt(abs(lat)),  $
                          count_bad)
          if (count_bad gt 0) then begin
              print, count_bad," Bad Point(s) found in dB/dt check"
              quality(loc_good(loc_bad)) = 1
          endif
      endif

      loc_missing = where(b(iComp,*) lt -70000,count)
      if (count gt 0) then begin
          quality(loc_missing) = -1
      endif

  endfor

  newb = fltarr(3, nPts, nDays)
  newq = fltarr(   nPts, nDays)

;  for iComp = 0,2 do begin
      for iDay = 0, nDays-1 do begin
          newq(      *,iDay) = quality(iDay*nPts:(iDay+1)*nPts-1)
      endfor
;  endfor

  quality = newq

end
