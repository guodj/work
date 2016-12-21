
pro quiet_day, data, time, lat, quality, quietdata, QuietDay, $
               UseQuietDay = UseQuietDay

  if (n_elements(UseQuietDay) eq 0) then UseQuietDay = 0

  nPts  = n_elements(data(0,*,0))
  nDays = n_elements(data(0,0,*))

  quietdata = data
  qd  = fltarr(3, nPts)
  qdm = fltarr(3, nPts)
  nd  = intarr(   nPts)

  if (UseQuietDay) then begin

      qdm = QuietDay

  endif else begin

      for iPt = 0, nPts-1 do begin

          loc = where(quality(iPt, *) eq 0, count)
          nd(iPt) = count
          if (count gt 0) then begin
              for iComp = 0, 2 do begin
                  qd(iComp, iPt) = median(data(iComp, iPt, loc))
              endfor
          endif

      endfor

      if (abs(lat) lt 60.0) then begin
          
          n8 = nPts/8 / 2

          for iPt = 0, nPts-1 do begin
              iStart = max([iPt-n8,0])
              iEnd   = min([iPt+n8,nPts-1])
              l = where(nd(iStart:iEnd) gt 0, count)
              if (count gt 0) then $
                for iComp=0, 2 do $
                qdm(iComp,iPt) = mean(qd(iComp,iStart+l))
          endfor

      endif else begin

          l = where(nd gt 0, count)
          if (count gt 0) then $
            for iComp=0, 2 do $
            qdm(iComp,l) = median(qd(iComp,l))

      endelse

      for iComp=0, 2 do begin
          lnz = where(qdm(iComp,*) ne 0, cnz)
          lz  = where(qdm(iComp,*) eq 0, cz)
          if (cz gt 0 and cnz gt 0) then qdm(iComp,lz) = mean(qdm(iComp,lnz))
      endfor

      QuietDay = qdm

  endelse

  for iPt = 0, nPts-1 do begin

      loc = where(quality(iPt, *) eq 0, count)
      nd(iPt) = count

      if (count gt 0) then begin
          for iComp = 0, 2 do begin
              quietdata(iComp, iPt, loc) = $
                data(iComp, iPt, loc) - qdm(iComp,iPt)
          endfor
      endif

      loc = where(quality(iPt, *) ne 0, count)
      if (count gt 0) then begin
          for iComp = 0, 2 do quietdata(iComp, iPt, loc) = -1.0e32
      endif

  endfor

end
