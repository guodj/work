
pro interpolate_ovation, ValueIn, LatsIn, MltsIn, $
                         LatsOut, MltsOut, ValueOut

  nLats = n_elements(LatsOut)
  nMlts = n_elements(MltsOut)

  ValueOut = fltarr(nMlts, nLats)

  for iMlt=0,nMlts-1 do begin

     ; Mlts are bounded from 0 to 24, so this should never fail...
     l = where(MltsIn ge MltsOut(iMlt))
     iM = l(0)
     if (iM eq 0) then begin
        iM  = 1
        xM = 1.0
     endif else begin
        xM = (MltsOut(iMlt) - MltsIn(iM)) / (MltsIn(iM) - MltsIn(iM-1))
     endelse

;     print, "mlt : ",iMlt, iM, xM, MltsOut(iMlt), MltsIn(iM)

     for iLat = 0,nLats-1 do begin

        l = where(LatsIn ge LatsOut(iLat), c)

        if (c gt 0 and l(0) gt 0) then begin
           iL = l(0)
           xL = (LatsOut(iLat) - LatsIn(iL)) / (LatsIn(iL) - LatsIn(iL-1))
        endif else begin
           iL  = 1
           xL = 1.0
        endelse
;        print, "lat : ",iLat, iL, xL, LatsOut(iLat), LatsIn(iL)

        ValueOut(iMlt,iLat) = $
           (1.0 - xM) * (1.0 - xL) * ValueIn(iM  ,iL  ) + $
           (      xM) * (1.0 - xL) * ValueIn(iM-1,iL  ) + $
           (1.0 - xM) * (      xL) * ValueIn(iM  ,iL-1) + $
           (      xM) * (      xL) * ValueIn(iM-1,iL-1)
           
     endfor

  endfor

end

