
pro imf_read, file, time, mag, vel, den, temp, nPts, notes

  notes = ''

  spawn, 'wc '+file,wc
  nLinesMax = long(wc(0))
  openr,1,file

  done = 0
  line = ''

  while (not done) do begin

      readf,1,line

      if eof(1) then done = 1  

      if (not done and strpos(line,'#START') eq -1) then notes = [notes,line]

      if (strpos(line,'#START') gt -1) then begin

          vel  = fltarr(3,nLinesMax)
          mag  = fltarr(3,nLinesMax)
          den  = fltarr(nLinesMax)
          temp = fltarr(nLinesMax)
          time = dblarr(nLinesMax)

          tmp = fltarr(15)
          itime = intarr(6)

          npts = 0L

          while (not eof(1)) do begin

              readf,1,tmp

              itime = tmp(0:5)
              c_a_to_r, itime, rtime
              time(npts)  = rtime
              mag(*,npts) = tmp(7:9)
              vel(*,npts) = tmp(10:12)
              den(npts)   = tmp(13)
              temp(npts)  = tmp(14)
              
              npts = npts + 1

          endwhile

          done = 1

      endif

  endwhile

  loc = where(vel(0,*) gt 0.0,count)
  if count gt 0 then vel(0,loc) = -1.0e32

  loc = where(den lt 0.0,count)
  if count gt 0 then den(loc) = -1.0e32

  loc = where(temp lt 0.0,count)
  if count gt 0 then temp(loc) = -1.0e32

  close,1

  if (n_elements(npts) eq 0) then begin
      print, "#START not found!!!"
      print, "Check imf file."
      stop
  endif

  time = time(0:npts-1)
  mag = mag(*,0:npts-1)
  vel = vel(*,0:npts-1)
  den = den(0:npts-1)
  temp = temp(0:npts-1)

end

