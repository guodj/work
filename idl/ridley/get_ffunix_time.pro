function get_ffunix_time, id, mrow

  common ffinfo, header

  point_lun, header.unit(id), long(mrow)*header.rl(id)

  if (eof(header.unit(id))) then begin
      t = -1.0e32
  endif else begin
      t = double(0.0)
      readu, header.unit(id), t
  endelse

  return, t

end
