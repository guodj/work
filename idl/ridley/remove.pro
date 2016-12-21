; procedure remove
;
;  Reads in one ascii file and outputs another ascii file which is identical
;  to the first minus all of the black spaces at the end of each line.
;

pro remove, filer$, filew$

  if n_elements(filer$) eq 0 and n_elements(filew$) eq 0 then begin

    filer$=''
    filew$=''
    read, 'enter file name:', filer$
    read, 'enter file to put into:', filew$

  endif

  openr,1,filer$
  openw,2,filew$

  dumb='' 

  while not eof(1) do begin

    readf,1,dumb

    printf,2, strtrim(dumb)

  endwhile

  close,1,2

  return

end
