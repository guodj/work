pro plot_data_new, datafile, date, time, dt, minlat, nodata=nodata, ut=ut

  if n_elements(nodata) eq 0 then nodata = 0 else nodata = 1

  scale_nt = 25.0
  scale_mv =  3.0

  !p.symsize=0.3

  mr = 90.0 - minlat

  if n_elements(ut) eq 0 then begin

    sdate = strmid(date,4,2)+'-'+strmid(date,0,3)+'-'+strmid(date,10,2)+' '
    sdate = sdate + strmid(time,0,5)

    c_s_to_a, itime, sdate
    c_a_to_r, itime, rtime

  endif else rtime = ut

  stime = rtime - dt
  etime = stime + 2.0*dt

  openr,11, datafile

  line = ''
  itime = intarr(6)

  done = 0

  while not done do begin

    readf, 11, line

    if eof(11) then begin
      print, 'EOF in datafile ',datafile
      stop
    endif

    if (strpos(line,'#TIME') gt -1) then begin
      readf, 11, itime
      c_a_to_r, itime, rtime_amie
      if (rtime_amie ge rtime) then done = 1
    endif

  endwhile

  if (rtime_amie gt rtime) then begin
    print, 'Mismatch in times in data file:'
    c_r_to_a, itime, rtime
    c_a_to_s, itime, stime
    print, '  Looking for time : ',stime
    c_r_to_a, itime, rtime_amie
    c_a_to_s, itime, stime
    print, '  Found time : ',stime
  endif

  done = 0

  ahn = 1
  mag = 2
  efield = 3
  fuv = 4

  data = fltarr(6)

  while (done eq 0) do begin

    type = 0

    while (type ge 0) do begin

      readf,11,line

      if eof(11) then type = -1

      if (strpos(line,'#TIME') gt -1) then type = -1
      if (strpos(line,'#AHN') gt -1) then type = ahn
      if (strpos(line,'#MAGNETOMETER') gt -1) then type = mag
      if (strpos(line,'#EFIELD') gt -1) then type = efield
      if (strpos(line,'#CONDUCTANCES') gt -1) then type = fuv

      if (type eq -1) then done = 1

      data = fltarr(6)
      if (type eq ahn) then psym = 2
      if (type eq mag) then psym = 4
      if (type eq efield) then begin
          psym = 5
          data = fltarr(4)
      endif
      if (type eq fuv) then begin
          psym = 1
          data = fltarr(7)
      endif

      if type gt 0 then begin

          if (type eq efield or type eq fuv) then readf,11,nFiles else nFiles=1

          for iFile = 0,nFiles-1 do begin

              if (type eq efield or type eq fuv) then readf,11,line

              readf, 11, npts

              for i = 0, npts-1 do begin 

                  readf, 11, lat, mlt, data
        
                  if (abs(lat) gt minlat) then begin

                      t = mlt*2.0*!pi/24.0 - !pi/2
                      r = 90.0 - abs(lat)

                      x0 = r*cos(t)
                      y0 = r*sin(t)

                      plots, [x0], [y0], psym = psym, symsize=0.3

                      if (type eq mag) then begin

                          h = data(0)
                          he = data(1)
                          e = data(2)
                          ee = data(3)

                          if (abs(h) lt 2000.0) and (abs(e) lt 2000.0) and $
                            abs(he) gt 0.0 and abs(ee) gt 0.0 then begin

                              theta = (2.0*!pi - t) mod (2.0*!pi)
                              d = (h^2.0 + e^2.0)^0.5

; want pure north (h) to be 0 degrees, and rotate clock wise is positive

                              phi = -1.0*asin(e/d)
                              if (h lt 0.0) then phi = !pi - phi

                              phi = !pi + (phi - theta)

; want to rotate counter clock wise 90 deg to the direction of convection

                              phi = phi + !pi/2.0

                              x1 = x0 + d*cos(phi)/scale_nt
                              y1 = y0 + d*sin(phi)/scale_nt

                              plots, [x0,x1], [y0,y1], thick = 3

                          endif

                      endif

                  endif

              endfor

          endfor

          type = 0

      endif

      if eof(11) then begin
        type = -1
        done = 1
      endif

    endwhile

  endwhile

  close,11

  return

end
