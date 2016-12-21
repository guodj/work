pro io_curve, x_final, y_final, in=in, out=out, screen = screen

  if n_elements(in) eq 0 then in = 0
  if n_elements(out) eq 0 then out= 0
  if n_elements(screen) eq 0 then screen = 0
  if out or in then begin

    if in then title = 'Enter File to Input'			$
    else title = 'Enter File to Output'

    if not screen then begin
      nf = 1
      mnf = 1
      get_files, filename, nf, mnf, title
      iofile = filename(0)
    endif else begin
      iofile = ''
      read, iofile
    endelse

    close, 1

  endif

  if in then begin

    openr,1,iofile

    n = 0
    readf,1,n

    x_final = fltarr(n)
    y_final = fltarr(n)
    x = 0.0
    y = 0.0

    for i=0,n-1 do begin

      readf,1, x, y

      x_final(i) = x
      y_final(i) = y

    endfor

    close, 1

  endif

  if out then begin

    openw,1,iofile

    n = n_elements(x_final)

    printf,1,n

    for i=0,n-1 do printf,1, x_final(i), y_final(i)

    close,1

  endif

  return

end

