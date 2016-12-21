
pro netcdfread, file

  id = ncdf_open(file)
  result = ncdf_inquire(id)
  help, result, /struct
  nvars = result.nvars
  natts  = result.ngatts

  for i=0,nvars-1 do begin
      r = ncdf_varinq(id,i)
      help, r, /struct
  endfor

  for i=0,natts-1 do begin
      name = ncdf_attname(id,i,/global)
;      att = ncdf_attinq(id, /global, name)
      ncdf_attget, id, /global, name, att
      print, name
      print, string(att)
  endfor

  ncdf_close, id

end
