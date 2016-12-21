pro write_tgcm_bin, fileout, data, vars, utdes, axdes, xax, yax, ntmax

  openw,11,fileout

  nv = n_elements(data(*,0,0))
  nx = n_elements(data(0,*,0))
  ny = n_elements(data(0,0,*))

  printf,11,nv,nx,ny
  for i=0, nv-1 do printf,11,vars(i)
  for i=0, nv-1 do printf,11,utdes(i)
  for i=0, nv-1 do for j=0,1 do printf,11,axdes(j,i)
  writeu,11,xax
  writeu,11,yax
  writeu,11,data

  close,11

  return

end

