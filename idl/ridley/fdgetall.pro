pro fdgetall, flatfile, time, data

  Setup_Var_Menu, 1, flatfile, timeint, vars, units, ncol, nrows, rowlen
  nr = nrows(0)
  nc = ncol(0)+1

  data = fltarr(nc,nr)
  time = dblarr(nr)
  dumb = fltarr(nc)
  dtime = double(0.0)

  openr, 11, flatfile+'.dat'

  i = long(0)

  for i=long(0),long(nr-1) do begin

    readu, 11, dtime
    readu, 11, dumb

    time(i) = dtime
    data(*,i) = dumb

  endfor

  close,11

  return

end

