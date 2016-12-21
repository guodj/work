
pro average_mags, data, time, lat, quality, nPtsMax, $
                  averagedata, averagetime, averagequality

  nPts  = n_elements(data(0,*,0))
  nDays = n_elements(data(0,0,*))

  averagedata = fltarr(3, nPtsMax, nDays) - 1.0e32
  averagetime = dblarr(   nPtsMax, nDays)

  nPtsDiv = nPts/nPtsMax

  averagequality = intarr(nPtsMax, nDays)

  for iDay = 0, nDays-1 do begin

      mt = mean(time(*,iDay))
      c_r_to_a, itime, mt
      itime(3:5) = 0
      c_a_to_r, itime, basetime
      dt = 24.0*3600.0 / nPtsMax

      print, "Averaging Day : ",iDay+1
      for iPt = 0, nPtsMax-1 do begin

          if (iPt eq 0) then t1 = basetime + dt * iPt $
          else t1 = basetime + dt * iPt - dt/2
          t2 = t1 + dt

          loc = where(time(*,iDay) ge t1 and $
                      time(*,iDay) lt t2 and $
                      quality(*, iDay) eq 0, count)
   
          if (count gt 0) then begin
              for iComp = 0, 2 do $
                averagedata(iComp, iPt, iDay) = $
                mean(data(iComp, loc, iDay))
              averagetime(iPt, iDay) = $
                mean(time(loc, iDay))
          endif else begin
              averagetime(iPt,iDay) = t1
              averagequality(iPt,iDay) = -1
          endelse
              
      endfor

  endfor


end

