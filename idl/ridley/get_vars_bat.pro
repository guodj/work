pro get_vars_bat, varname, varcount, numofvar

  done = 0
  nf = n_elements(varcount)

  while done eq 0 do begin

    que = ''
    print, 'By entering numbers, you will select or deselect var'
    print, 'Enter file number, column number (return to exit) :'
    read, que

    if strlen(que) eq 0 then done = 1 else begin

      loc = where(byte(que) lt 46 or byte(que) gt 57, count)
      if count gt 0 then begin
	find = where(byte(que) ge 46 and byte(que) le 57, num)
	if num gt 0 then begin
	  r = fix(que)-1
	  c = fix(strmid(que,loc(count-1)+1,strlen(que)-loc(count-1)-1)) - 2
	endif else begin
	  print, 'Enter numerical values, please'
	  r = nf + 1
	endelse
      endif else begin
	print, 'Enter 2 numbers (file,column) !'
	r = nf+1
      endelse

      if r le nf then begin

        if (c le varcount(r)) and (c ge 0) then begin

	  ev = varname(r,c)

	  if strmid(ev,14,10) eq '(selected)' then begin
	    numofvar = numofvar-1
	    ev = strmid(ev,0,14)
	    print, ev,' - Deselected'
	  endif else begin
	    numofvar = numofvar + 1
	    ev = strmid(ev,0,14) + '(selected)'
	    print, ev
	  endelse

	  varname(r,c) = ev

	endif

	if (c gt varcount(r)) or (c lt -1) then print, 'Column is out of range'
	if (c eq -1) then print, 'Time is not an option!!!'

      endif

    endelse

  endwhile
  print, 'You now have ',tostr(numofvar+1),' variables selected.'

  return

end  
