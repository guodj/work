pro get_station, statfile, xyzpos, stat, lat, lon, alpha, mln

  dum$=''
  xyzpos = intarr(1,3)

  close, 1
  openr, 1, statfile(0)
  k = 0

  while not eof(1) do begin

    readf, 1, dum$
    if xyzpos(0,0) eq 0 then begin
      stat = [strmid(dum$,0,3)]
      lon = [float(strmid(dum$,10,4))]
      lat = [float(strmid(dum$,15,5))]
      alpha = [float(strmid(dum$,21,3))]
      mln = [float(strmid(dum$,25,4))]
    endif else begin
      stat = [stat,strmid(dum$,0,3)]
      lon = [lon,float(strmid(dum$,10,4))]
      lat = [lat,float(strmid(dum$,15,5))]
      alpha = [alpha,float(strmid(dum$,21,3))]
      mln = [mln,float(strmid(dum$,25,4))]
      xyzpos = [xyzpos,intarr(1,3)]
    endelse
    xyzpos(k,0) = fix(strmid(dum$,7,2))
    xyzpos(k,1) = fix(strmid(dum$,4,2))
    xyzpos(k,2) = 2*xyzpos(k,0) - xyzpos(k,1)
    k = k + 1

  endwhile

  close, 1

  RETURN

END

