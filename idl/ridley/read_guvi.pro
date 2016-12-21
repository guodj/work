
pro read_guvi, file, time, data, latitude, longitude

  id = ncdf_open(file)

  ncdf_varget, id, 'PiercePointLatitude', latitude
  ncdf_varget, id, 'PiercePointLongitude', longitude
  ncdf_varget, id, 'RadianceData', data
  ncdf_varget, id, 'Time', utime
  ncdf_varget, id, 'DOY', doy

  ncdf_attget, id, /global, "Starting_Time", stime
  stime = string(stime)

  ncdf_close, id

  npts = n_elements(utime)
  time = dblarr(npts)

  year = fix(strmid(stime,0,4))

  for i=0,npts-1 do begin
      if (i gt 0) then if (doy(i) lt doy(i-1)) then year = year + 1
      itime = [year, 1, doy(i), 0,0,0]
      c_a_to_r, itime, rtime
      time(i) = rtime + double(utime(i))/1000.0
  endfor

end
