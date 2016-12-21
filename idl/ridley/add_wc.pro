pro add_wc, filelist

  if n_elements(filelist) eq 0 then begin

    filelist = strarr(1)

    while strlen(filelist(0)) eq 0 do begin

      filein = ''
      print, 'Enter file or group of files (*.ave, ect) to add word count :'
      read, filein
      filelist = findfile(filein)

    endwhile

  endif

  n = n_elements(filelist)

  for i=0,n-1 do print, 'File # ',strcompress(string(i)), ' - ',filelist(i)

  que = ''
  read, 'Is this list correct (y/n) ? ', que

  if (strmid(que,0,1) eq 'y') or (strmid(que,0,1) eq 'Y') then begin

    for i=0,n-1 do begin

      spawn, 'wc '+filelist(i),numlines
      nl = fix(numlines) - 2
      snl = '     '+strcompress(string(nl))
      snl = strmid(snl,strlen(snl)-5,5)

      close,1
      openu,1,filelist(i)

      line = ''
      readf,1,line
      line = snl+strmid(line, 4,strlen(line)-5)

      printf,1, line

      close, 1

    endfor

  endif

  return

end
