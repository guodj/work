
pro readdmsp, dmspfile, date, time, maxran, dmspdata, 		$
	interpolate=interpolate

  close,1
  openr,1,dmspfile

  line = ''
  readf,1,line
  readf,1,line

  stime = strmid(date,4,2)+'-'+mkupper(strmid(date,0,3))+'-'+strmid(date,10,2)
  stime = stime+' '+strmid(time,0,5)
  c_s_to_a,itime,stime
  c_a_to_r,itime,amietime
  dt = 30.0*60.0

  on = 0
  done = 0
  pot = 0.0
  dmspdata = fltarr(6,1000)
  n = 0
  dold = 1.0

  while (not(eof(1)) and not(done)) do begin
    readf,1,line
    alat = float(strmid(line,19,6))
    if alat gt 90.0-maxran then begin
      year = fix(strmid(line,0,4)) 
      mo = fix(strmid(line,5,2))
      da = fix(strmid(line,8,2))
      hr = fix(strmid(line,11,2))
      mn = fix(strmid(line,14,2))
      sc = fix(strmid(line,17,2))
      itime = [year,mo,da,hr,mn,sc] 
      c_a_to_r,itime,dmsptime
      if (abs(amietime-dmsptime) lt dt) then begin
        alon = float(strmid(line,26,6))
        mlt = float(strmid(line,33,5))
        kloc = float(strmid(line,39,4))
        ang1n = float(strmid(line,43,7))
        ealong = float(strmid(line,50,7))
        ste = float(strmid(line,57,7))
        bgauss = float(strmid(line,77,8))
        glat = float(strmid(line,107,6))
        glon = float(strmid(line,114,6))

        the = mlt*!pi/12.0 - !pi/2.0
        x = (90.0-alat)*cos(the)
        y = (90.0-alat)*sin(the)
	if (on eq 1) then begin
	  d = distance(glat, glon, lato, lono)
	  if (d gt 500.0) then d = dold
	  pot = pot - 				$
		ealong*d/1000.0
	  dmspdata(0,n) = alat
	  dmspdata(1,n) = mlt
	  dmspdata(2,n) = pot
	  dmspdata(3,n) = x
	  dmspdata(4,n) = y
          dmspdata(5,n) = float(hr) + float(mn)/60.0 + float(sc)/3600.0
	  n = n+1
        endif
	on = 1
        lato = glat
        lono = glon
	if (n_elements(d) gt 0) then dold = d
      endif else if (on eq 1) then done = 1
    endif
  endwhile

  close,1

  dmspdata = dmspdata(0:5,0:n-1)

  if n_elements(interpolate) gt 0 then begin
    endp = dmspdata(2,n-1)
    for i=0,n-1 do begin
      dmspdata(2,i) = dmspdata(2,i) - float(i)*endp/float(n)
    endfor
  endif

end
