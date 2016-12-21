pro convert_mhd_ascii_to_save, input

  list = findfile(input)
  nfiles = n_elements(list)

  for l=0,nfiles-1 do begin

    file = list(l)

    print, 'Reading from file : ',file
    openr,1,file

    line = ''
    for i=0,8 do readf,1,line
    nlats = fix(strmid(line,strpos(line,'I=')+2,4))
    nlons = fix(strmid(line,strpos(line,'J=')+2,5))

    nvars = 18

    if l eq 0 then begin
      tmp = fltarr(nvars)
      data = fltarr(2,nvars,nlons,nlats)
    endif

    for n=0,1 do begin
      for j=0,nlons-1 do for i=0,nlats-1 do begin
        readf,1,tmp
        data(n,*,j,i) = tmp
      endfor
      if n eq 0 then readf,1,line
    endfor

    close,1

    print,'Saving to file ',file+'.save'
    save, data,nvars,nlons,nlats, file=file+'.save'

  endfor

  return

end

