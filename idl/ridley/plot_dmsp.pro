
filelist = findfile('a139????.*')

nfiles = n_elements(filelist)

ppp = 9
space = 0.01
pos_space, ppp, space, sizes
total_plots = sizes.nby*sizes.nbx

n_start = 0
n_end = nfiles-1

for n=n_start, n_end do begin

  pn = n mod ppp
  if pn eq 0 then plotdumb

  openr,1,filelist(n)

  tmp=fltarr(19)

  readf,1,tmp

  satid = tmp(0)
  itime = fltarr(6)
  itime(0) = tmp(1)
  itime(1) = 1
  itime(2) = tmp(2)
  itime(3) = 0
  itime(4) = 0
  itime(5) = 0

  c_a_to_r, itime, basetime

  nlines = fix(tmp(16))

  data = fltarr(13, nlines)

  line = ''
  for i=0,nlines-1 do begin
    readf,1,line
    data(0,i) = float(strmid(line,0,7))
    data(1,i) = float(strmid(line,7,7))
    data(2,i) = float(strmid(line,14,7))
    data(3,i) = float(strmid(line,21,7))
    data(4,i) = float(strmid(line,28,7))
    data(5,i) = float(strmid(line,35,3))
    data(6,i) = float(strmid(line,38,8))
    data(7,i) = float(strmid(line,46,7))
    data(8,i) = float(strmid(line,53,7))
    data(9,i) = float(strmid(line,60,7))
    data(10,i) = float(strmid(line,67,7))
    data(11,i) = float(strmid(line,74,8))
    data(12,i) = float(strmid(line,81,8))
  endfor

  close,1

  time = basetime + reform(data(0,*))

  vy = reform(data(1,*))
  vz = reform(data(2,*))

  pot = reform(data(6,*))
  lat = reform(data(8,*))
  mlt = reform(data(9,*))/24.0 * 2.0 * !pi - !pi/2.0

  get_position, ppp, space, sizes, pn, pos

  mr = 40.0

  plot, [-mr,mr], [-mr,mr], xstyle = 5, ystyle = 5, pos = pos, $
	/nodata, /noerase

  if mean(lat) lt 0.0 then begin
    title = 'South'
    lat = 90.0 + lat
    yfac = -1.0
    xfac = 1.0
  endif else begin
    title = 'North'
    lat = 90.0 - lat
    yfac = 1.0
    xfac = -1.0
  endelse

  x = lat * cos(mlt)
  y = lat * sin(mlt)

  dx = x(1:nlines-1) - x(0:nlines-2)
  dy = y(1:nlines-1) - y(0:nlines-2)

  dxnew = dx
  dynew = dy

  nsmooth = 10

  for i = nsmooth/2, nlines - nsmooth/2 -2 do begin

    dxnew(i) = mean(dx(i-nsmooth/2:i+nsmooth/2-1))
    dynew(i) = mean(dy(i-nsmooth/2:i+nsmooth/2-1))

  endfor

  dx = dxnew
  dy = dynew

  for i=0, nlines - 2 do begin

    if (lat(i) lt mr and vy(i) lt 100.0) then begin

      dr = sqrt(dx(i)^2 + dy(i)^2)

      vel_y = 10.0 * vy(i) * dx(i) / dr * yfac
      vel_x = 10.0 * vy(i) * dy(i) / dr * xfac

      oplot, [x(i), x(i) + vel_x], [y(i), y(i) + vel_y]

    endif

  endfor

  no00 = 1
  no06 = 1
  no12 = 1
  no18 = 1

  if (n_end-n lt sizes.nbx) or 						$
     (pn+1 gt total_plots-sizes.nbx) then no00 = 0
  if (pn mod sizes.nbx eq sizes.nbx-1) or 			$
     (n_end-n eq 0) then no06 = 0
  if (pn lt sizes.nbx) then no12 = 0
  if (pn mod sizes.nbx eq 0) or (n eq n_start) then no18 = 0

;--------------------------------------------------------------
; Draw the MLT grid
;--------------------------------------------------------------

  plotmlt, mr, no00 = no00, no06 = no06, no12 = no12, no18 = no18

  xyouts, pos(2), pos(3), title, /norm, align=1.0

  loc = where(lat eq min(lat))
  ut = time(loc(0))
  c_r_to_a, itime, ut
  c_a_to_s, itime, stime

  date = strmid(stime, 10,5)

  print, stime

  xyouts, pos(0), pos(3), date, /norm

  c_r_to_a, itime, time(loc(0))-3600.0

  dir = '/g/ridley/BATSRUS_o2k/run.'+tostr(itime(0))+ $
	chopr('0'+tostr(itime(1)),2)+ $
	chopr('0'+tostr(itime(2)),2)+'.'+ $
	chopr('0'+tostr(itime(3)),2)+ $
	chopr('0'+tostr(itime(4)),2)+ $
	chopr('0'+tostr(itime(5)),2)+'/ionosphere'

  print, dir
  spawn, 'cp '+filelist(n)+' '+dir

endfor

end
