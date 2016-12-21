
pro calc_efield, potential, mlt, lat, efielde, efieldn

  nLats = n_elements(potential(0,*))
  nMlts = n_elements(potential(*,0))

  EFieldN = fltarr(nMlts, nLats)
  EFieldE = fltarr(nMlts, nLats)

  if (max(abs(potential)) lt 1000.0) then fac = 1000.0 else fac = 1.0

  for iMlt = 1, nMlts-2 do for iLat = 1, nLats-2 do begin

      EFieldN(iMlt, iLat) = - fac * $
         (Potential(iMlt,iLat+1) - Potential(iMlt,iLat-1)) / $
        ((Lat(iLat+1) - Lat(iLat-1))*!pi/180.0 * 6372000.0)

      EFieldE(iMlt, iLat) = - fac * $
         (Potential(iMlt+1,iLat) - Potential(iMlt-1,iLat)) / $
        ((MLT(iMlt+1) - MLT(iMlt-1))*!pi/12.0 * 6372000.0 * $
         cos(Lat(iLat)*!pi/180.0))

  endfor

  EFieldE(0,*) = (EFieldE(1,*)+EFieldE(nMlts-2,*))/2.0
  EFieldE(nMlts-1,*) = EFieldE(0,*)

  EFieldN(0,*) = EFieldN(1,*)
  EFieldN(nLats-1,*) = EFieldN(nLats-2,*)

end

