function get_file_batch,max

  infile = strarr(max)
  print, 'You can enter up to ',tostr(max),' file names. Please press'
  print, 'return when you are done entering file names.'

  numoffil=0

  done = 0
  while (done eq 0) do begin

    que$ = ''
    print, 'Enter file name ',strcompress(string(numoffil+1),/remove),' :'
    read, que$
    if (strmid(que$,0,1) eq ' ') or (strlen(que$) eq 0) then done = 1 	$
    else begin
      infile(numoffil) = que$
      numoffil = numoffil+1
    endelse
    if numoffil ge max then done = 1

  endwhile

  info = {files : infile, nf : numoffil}

  return, info

end
