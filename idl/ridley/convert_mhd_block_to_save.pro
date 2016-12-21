
pro convert_mhd_block_to_save, filelist

  list = findfile(filelist)
  n = n_elements(list)
  if (n eq 1) and (strlen(list(0)) eq 0) then begin
    print, 'Sorry, not files found'
  endif else begin

    line = ''
    ncols = 23

    for i=0,n-1 do begin

      file = list(i)
      print, 'Working on file : ',file

      zipped = 0
      s = strpos(file, '.gz')
      if s gt -1 then begin
        zipped = 1
        print, 'Unzipping file.'
        spawn, 'gunzip '+file
	file = strmid(file, 0, s)
      endif

      n_equals = -1

      openr,1,file
      n = 0
      while (n_equals lt 0) do begin
        readf,1,line
        npos = strpos(line,"N=")
        if (npos gt -1) then begin
          nrows = float(strmid(line,npos+2,10))
          n_equals = n
        endif
	n = n + 1
      endwhile
      close,1

      tmp = fltarr(nrows)
      data = fltarr(ncols, nrows)

      openr,1,file
      for n=0,n_equals do begin
	readf,1,line
      endfor

      print, 'Reading Data.'
      for n=0,ncols-1 do begin
        readf,1,tmp
	data(n,*) = tmp
      endfor

      close,1

      save, data, file=file+'.save'

      if (zipped) then begin
	print, 'Rezipping the file.'
        spawn, 'gzip '+file
      endif

    endfor

  endelse

  return

end