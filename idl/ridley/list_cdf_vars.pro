
pro list_cdf_vars, file, vars

  id = cdf_open(file)
  r = cdf_inquire(id)
  nVars = r.nVars
  if (nVars eq 0) then begin
     nVars = r.nzVars
     zvars = 1
  endif else zvars = 0

  vars = strarr(nVars)
  for i=0,nVars-1 do begin

     if (zvars) then var = cdf_varinq(id,i,/zvariable)
     if (not zvars) then var = cdf_varinq(id,i)
     
     vars(i) = var.name

  endfor

end

