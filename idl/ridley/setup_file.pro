pro setup_file, fbase, title, vars, n

  nele = n_elements(vars(n,*))

  fbase = widget_base(title=title,/row)

  dum = widget_button(fbase, uvalue='DONE', value = 'Done')

  bbase = widget_base(fbase, /row)

  nbase = lonarr(5)
  div = (nele-1)/5 + 1
  if div le 0 then div = 1

  for i=0,4 do nbase(i) = widget_base(bbase, /column)

  for i=0,nele-1 do begin

    if strlen(vars(n,i)) gt 0 then begin
      if strmid(vars(n,i),14,10) ne '(selected)' then		$
	vars(n,i) = strmid(vars(n,i),0,14)+'                '
      dummy = widget_button(nbase(i/div), value = vars(n,i), 	$
	uvalue=vars(n,i))
    endif

  endfor

  return

end
