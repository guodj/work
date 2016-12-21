pro printps, directory

  if n_elements(directory) eq 0 then begin

    directory = ''
    print, 'enter directory :'
    read, directory

  endif

  psfiles = findfile(directory+'*.ps')

  if strlen(psfiles(0)) eq 0 then begin

    print, 'There are no ps files in this directory!!'

  endif else begin

    for i=0,n_elements(psfiles)-1 do begin

      spawn, 'lpp '+psfiles(i)
      wait, 5.0

    endfor

  endelse

  return

end
