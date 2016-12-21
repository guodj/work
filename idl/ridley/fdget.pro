pro fdget, id, mrow, uttime, data

  common ffinfo, header

  point_lun, header.unit(id), long(mrow)*header.rl(id)

  uttime = double(0.0)
  nvar = (header.rl(id)-8)/4
  data = fltarr(nvar)

  readu, header.unit(id), uttime

  dumb = 0.0
  for i=0,nvar-1 do begin
    readu, header.unit(id), dumb
    data(i) = dumb
  endfor

  return

end
