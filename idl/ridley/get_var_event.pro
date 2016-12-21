pro get_var_event, ev, screen=screen

  common titles, var, nvar

  if n_elements(screen) eq 0 then 				$
    widget_control, ev.id, get_uvalue = value			$
  else								$
    value = ev

  if (n_elements(value) eq 0) then value=''

  if value eq 'DONE' then widget_control, /destroy, ev.top		$
  else begin

    for i=0,n_elements(var)-1 do if value eq var(i) then n=i
    if n_elements(n) gt 0 then begin

      if strmid(var(n),14,10) eq '(selected)' then begin
	nvar = nvar - 1
	var(n) = strmid(var(n),0,14)+'                '
      endif else begin
	nvar = nvar + 1
	var(n) = strmid(var(n),0,14)+'(selected)'
      endelse

      if n_elements(screen) eq 0 then 				$
        widget_control, ev.id, set_value=var(n), set_uvalue=var(n)

    endif

  endelse

  return

end 

