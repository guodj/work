function get_cdf_fill, id, var, zvar

  if zvar eq 0 then begin
    cdf_attget, id, 'FILLVAL', var, value
  endif else begin
    cdf_attget, id, 'FILLVAL', var, value, /zvariable
  endelse

  return, value(0)

end

