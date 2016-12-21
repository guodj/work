pro change_amie_resolution, amie_data, in_lats, in_mlts, out_lats, out_mlts

  nt = n_elements(amie_data(*,0,0))
  nl = n_elements(out_lats)
  nm = n_elements(out_mlts)

  new_lat_interp = (90.0-out_lats)/(in_lats(0)-in_lats(1))
  new_mlt_interp = out_mlts/(in_mlts(1)-in_mlts(0))

  lat2d = fltarr(nm,nl)
  lon2d = fltarr(nm,nl)

  for i=0,nl-1 do lon2d(*,i) = new_mlt_interp
  for i=0,nm-1 do lat2d(i,*) = new_lat_interp

  new_amie_data = fltarr(nt,nm,nl)

  for i=0,nt-1 do begin

    new_amie_data(i,*,*) = interpolate(reform(amie_data(i,*,*)),lon2d,lat2d)

  endfor

  amie_data = new_amie_data

  return

end

