
pro calculate_on2, oIn, n2In, alts, on2

  nAlts = n_elements(alts)
  nLons = n_elements(oIn(*,0,0))
  nLats = n_elements(oIn(0,*,0))

  n2  = fltarr(nLons,nLats)
  o   = fltarr(nLons,nLats)
  on2 = fltarr(nLons,nLats)

  MaxValN2 = 1.0e21

  for iLon = 0, nLons-1 do begin
     for iLat = 0, nLats-1 do begin

        iAlt = nAlts-1
        Done = 0
        if (max(n2In(iLon,iLat,*)) eq 0.0) then Done = 1
        while (Done eq 0) do begin
           dAlt = (Alts(iAlt)-Alts(iAlt-1))*1000.0
           n2Mid = (n2In(iLon,iLat,iAlt) + $
                    n2In(iLon,iLat,iAlt-1)) /2.0
           oMid  = (oIn(iLon,iLat,iAlt) + $
                    oIn(iLon,iLat,iAlt-1)) /2.0

           if (n2(iLon,iLat) + n2Mid*dAlt lt MaxValN2) then begin
              n2(iLon,iLat) = n2(iLon,iLat) + n2Mid*dAlt
              o(iLon,iLat)  =  o(iLon,iLat) +  oMid*dAlt
              iAlt = iAlt - 1
           endif else begin
              dAlt = (MaxValN2 - n2(iLon,iLat)) / n2Mid
              n2(iLon,iLat) = n2(iLon,iLat) + n2Mid*dAlt
              o(iLon,iLat)  =  o(iLon,iLat) +  oMid*dAlt
              Done = 1
           endelse
        endwhile

     endfor
  endfor

  loc = where(n2 gt 0.0,count)
  if (count gt 0) then on2(loc) = o(loc)/n2(loc)

end
