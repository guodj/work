
pro read_amie_nc, amie_file, data, lats, mlts, time

  id = ncdf_open("apr18_2002_nh_pot.nc")

  ncdf_varget, id, 'lat', lats
  ncdf_varget, id, 'lon', mlts
  mlts = mlts*24.0/360.0
  ncdf_varget, id, 'ut', time
  time = time*3600.0
  ncdf_varget, id, 'pot', pot

  nm = n_elements(pot(*,0,0))
  nl = n_elements(pot(0,*,0))
  nt = n_elements(pot(0,0,*))

  data = fltarr(nt,1,nm,nl)
  for i=0,nt-1 do data(i,0,*,*) = pot(*,*,i)/1000.0

end
