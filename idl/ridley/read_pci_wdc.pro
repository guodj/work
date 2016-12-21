
pro read_pci_wdc, filein, pci, time

  spawn, 'wc '+filein,output
  nlines = fix(output)

  openr,1,filein

  iline = 0

  line = ''

  while (not eof(1) and iline lt nlines) do begin

    readf,1,line

    if (iline eq 0) then begin
      n = strlen(line)
      npts = (n - 24) / 6
      pci = fltarr(npts*nlines)
      time = dblarr(npts*nlines)
    endif

    iyear  = fix(strmid(line,3,2))
    imonth = fix(strmid(line,5,2))
    iday   = fix(strmid(line,8,2))

    itime = [iyear, imonth, iday, 0, 0, 0]
    c_a_to_r, itime, starttime

    for i=0,npts-1 do begin
      n = iline*npts + i
      pci(n) = float(strmid(line,24+i*6,6))
      time(n) = starttime + double(i) * 15.0 * 60.0
    endfor

    iline = iline + 1

  endwhile

  close,1

  return

end