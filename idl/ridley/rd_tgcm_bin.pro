pro rd_tgcm_bin, filein, data, vars, utdes, axdes, xax, yax

  openr,11,filein

  nv = 0
  nx = 0
  ny = 0

  readf,11,nv,nx,ny

  vars  = strarr(nv)
  utdes = strarr(nv)
  axdes = strarr(2,nv)
  xax   = fltarr(nx)
  yax   = fltarr(ny)
  data = fltarr(nv,nx,ny)

  tmp = ''
  for i=0, nv-1 do begin
    readf,11,tmp
    vars(i) = tmp
  endfor
  for i=0, nv-1 do begin
    readf,11,tmp
    utdes(i)=tmp
  endfor
  for i=0, nv-1 do for j=0,1 do begin
    readf,11,tmp
    axdes(j,i)=tmp
  endfor

  readu,11,xax
  readu,11,yax
  readu,11,data

  close,11

  return

end

