
pro display_cdf_variables, filename

  id = cdf_open(filename)
  re = cdf_inquire(id)
  nvars = re.nvars  
  nrows = re.maxrec

  if (nrows lt 1) then begin
    cdf_control, id, variable=0, /zvariable,get_var_info=v
    nrows = v.maxrec
  endif

  print, nrows
  if nVars eq 0 then begin
      nVars = re.nzvars
      print, "Z Vars"
      for i=0,nvars-1 do begin
          result = cdf_varinq(id, i, /zvar)
          print, i,'. ',result.name, result.dim
;          help, result,/struct
      endfor
  endif else begin

      for i=0,nvars-1 do begin

          result = cdf_varinq(id, i)
          print, i,'. ',result.name

      endfor

  endelse

  cdf_close,id

end
