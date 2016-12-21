
setdevice,'idl.ps','p'

filelist = findfile('a???????.*')

nfiles = n_elements(filelist)

ppp = 12
space = 0.01
pos_space, ppp, space, sizes, nx = 3
total_plots = sizes.nby*sizes.nbx

n_start = 0
n_end = nfiles-1

rms_save = fltarr(2,n_end-n_start+1)
time_save = dblarr(n_end-n_start+1)
north_save = intarr(n_end-n_start+1)
f13_save = intarr(n_end-n_start+1)
cpcp_save = fltarr(2,n_end-n_start+1)

for n=n_start, n_end do begin

  pn = (n*3) mod ppp
  if pn eq 0 then plotdumb

  openr,1,filelist(n)

  tmp=fltarr(19)

  readf,1,tmp

  if tmp(0) eq 13 then f13_save(n-n_start) = 1

  satid = 'F'+tostr(fix(tmp(0)))
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

  underplot = 0

  if (not eof(1)) then begin

    readf,1,nx,ny
    mhd_pot_input = fltarr(nx,ny)
    xpos_input = fltarr(nx,ny)
    ypos_input = fltarr(nx,ny)
    readf,1,mhd_pot_input
    readf,1,xpos_input
    readf,1,ypos_input
    underplot = 1

  endif

  close,1

  time = basetime + reform(data(0,*))

  savetime = reform(data(0,*))

  mhd  = reform(data(1,*))
  dmsp = reform(data(2,*))

  pot  = reform(data(6,*))
  potd = reform(data(3,*))
  lat  = reform(data(8,*))
  mlt  = reform(data(9,*))/24.0 * 2.0 * !pi - !pi/2.0

  if mean(lat) lt 0.0 then satid = ' S '+satid else begin
    satid = ' N '+satid 
    north_save(n-n_start) = 1
  endelse

  get_position, ppp, space, sizes, pn, pos
  pos(0) = pos(0) - space*2
  pos(2) = pos(2) - space*2

  mr = 40.0

  if underplot then begin

    r = sqrt(xpos_input^2+ypos_input^2)
    loc = where(r lt mr)

    mini = min(mhd_pot_input(loc))
    maxi = max(mhd_pot_input(loc))

    mini = -100.0
    maxi = 100.0

    if (maxi eq mini) then maxi = mini + 1.0

    maxi = max([abs(mini),maxi])
    mini = -maxi
    levels = fltarr(16)
    levels(0:7) = -maxi*(8-findgen(8))/8.0
    levels(8:15) = maxi*(findgen(8)+1)/8.0
    readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    c_colors = (ncolors-1)*findgen(30)/29.0 + 1

    contour, mhd_pot_input(loc), xpos_input(loc), ypos_input(loc), /irr,    $
        /follow, nlevels=30, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,/cell_fill

    contour, mhd_pot_input(loc), xpos_input(loc), ypos_input(loc), /irr, $
        /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

    mini_tmp = min(pot)
    maxi_tmp = max(pot)
    cpcp_m = maxi_tmp - mini_tmp

    locdmsp = where(abs(potd) lt 500.0)

    mini_tmp = min(potd(locdmsp))
    maxi_tmp = max(potd(locdmsp))
    cpcp_d = maxi_tmp - mini_tmp

    mhd_s  = string(cpcp_m,format="(f5.1)")+" MHD"
    dmsp_s = string(cpcp_d,format="(f5.1)")+" DMSP"

    charsize = 0.85

    xyouts, pos(2)+3*space,pos(1)+space*2.4,'CPCP', /norm, $
	align=1.0, charsize = charsize
    xyouts, pos(2)+3*space,pos(1)+space*1.2,mhd_s, /norm, $
	align=1.0, charsize = charsize
    xyouts, pos(2)+3*space,pos(1),          dmsp_s, /norm, $
	align=1.0, charsize = charsize

  endif else begin

    plot, [-mr,mr], [-mr,mr], xstyle = 5, ystyle = 5, pos = pos, $
	/nodata, /noerase

  endelse

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

  istart = 0

  for i=0, nlines - 2 do begin

    if (lat(i) lt mr and dmsp(i) lt 100.0) then begin

      dr = sqrt(dx(i)^2 + dy(i)^2)

      vel_y = 10.0 * dmsp(i) * dx(i) / dr * yfac
      vel_x = 10.0 * dmsp(i) * dy(i) / dr * xfac

      oplot, [x(i), x(i) + vel_x], [y(i), y(i) + vel_y]

      if istart eq 0 then begin
        if (x(i) lt 0.0) then xyouts, x(i)-1.0, y(i), 'S', align=1.0
        if (x(i) gt 0.0) then xyouts, x(i)+1.0, y(i), 'S'
	istart = 1
      endif


    endif

  endfor

  no00 = 1
  no06 = 1
  no12 = 1
  no18 = 0

  if pn eq 0 then no12 = 0
  if pn eq ppp - 3 or n eq n_end then no00 = 0

;--------------------------------------------------------------
; Draw the MLT grid
;--------------------------------------------------------------

  plotmlt, mr, no00 = no00, no06 = no06, no12 = no12, no18 = no18

  p1 = pos(0) + (pos(2) - pos(0))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95
  p2 = pos(3) - (pos(3) - pos(1))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95

  loc = where(lat eq min(lat))
  ut = time(loc(0))
  c_r_to_a, itime, ut
  c_a_to_s, itime, stime

  time_save(n) = ut

  outtime = tostr(itime(0)) + $
	chopr('0'+tostr(itime(1)),2) + $
	chopr('0'+tostr(itime(2)),2) + '.' + $
	chopr('0'+tostr(itime(3)),2) + $
	chopr('0'+tostr(itime(4)),2) + $
	chopr('0'+tostr(itime(5)),2)

  xyouts, p1, p2, $
          strmid(stime,10,2)+strmid(stime,13,5)+' UT', $
          /norm, alignment = 0.5, charsize = 0.8, orient = 45.0

  p1 = pos(0) + (pos(2) - pos(0))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.65
  p2 = pos(3) - (pos(3) - pos(1))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.65

  xyouts, p1, p2, $
          strmid(stime,0,9), $
          /norm, alignment = 0.5, charsize = 0.8, orient = 45.0

  xyouts, pos(0)-space*2.0, (pos(1)+pos(3))/2.0, title, $
	/norm, align=0.5, orient = 90.0

  xyouts, pos(0)-space*2.0, pos(3), strmid(satid,3,3), $
	/norm, align=0.0

  get_position, ppp, space, sizes, pn+1, pos1
  get_position, ppp, space, sizes, pn+2, pos

  pos(0) = pos1(0)

  pos(0) = pos(0) + space*3
  pos(2) = pos(2) + space*2

  pos(1) = pos(1) + space*2

  loc = where(lat lt mr,count)

  if (count gt 0) then begin

    stime = min(time(loc))
    etime = max(time(loc))

    time_axis, stime, etime, s_time_range, e_time_range, 	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

    plottime = time(loc) - stime
    plotdmsp = dmsp(loc)
    plotmhd  = mhd(loc)
    potdd    = potd(loc)

    loc = where(abs(plotdmsp) lt 10.0 and potdd ne 0.0,count)
    
    yrange = mm([plotdmsp(loc),plotmhd(loc),-2.0,2.0])

    if pn eq 0 then title = 'Velocity (km/s)' else title = ''

    plot, plottime(loc), plotdmsp(loc), yrange = yrange, $
        xrange = [s_time_range, e_time_range], xtickname = xtickname, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        pos = pos, /noerase, max_val = 10.0, title = title

    oplot, plottime, plotmhd, linestyle = 2
    oplot, plottime, plottime*0.0, linestyle = 1

    xpos = (pos(2) - pos(0))*0.02 + pos(0)
    ypos = (pos(3) - pos(1))*0.10 + pos(1)
    if count gt 0 then begin
      pot_rms = sqrt(mean((plotdmsp(loc) - plotmhd(loc))^2.0))
      tot_rms = sqrt(mean((plotdmsp(loc))^2.0))
    endif

    rms_save(0,n) = pot_rms
    rms_save(1,n) = tot_rms

    cpcp_save(0,n) = max(mhd_pot_input) - min(mhd_pot_input)
    cpcp_save(1,n) = max(potd(loc)) - min(potd(loc))

    xyouts, xpos, ypos, 'RMS : '+tostrf(pot_rms), $
        charsize = 0.8, /norm
    ypos = (pos(3) - pos(1))*0.04 + pos(1)
    xyouts, xpos, ypos, 'RMS : '+tostrf(tot_rms)+' (DMSP Only)', $
        charsize = 0.8, /norm

    print, format = '(a,a,2f9.4)',outtime,satid,tot_rms,pot_rms

    xyouts, pos(2)+space, pos(3), filelist(n), orient = -90.0, /norm,$
	charsize = 0.6

  endif

;  openw,3,filelist(n)+'.final'

  npts = n_elements(savetime)

  istart = 0
  iend   = 0

  for i=0,npts-1 do begin

    if (potd(i) ne 0.0) then istart = istart + 1

    if (potd(i) eq 0.0) then begin
      iend = 1
      for j=i,npts-1 do if (potd(j) ne 0.0) then iend = 0
    endif

;    if (istart gt 0 and iend eq 0) then printf,3,savetime(i),mhd(i)

  endfor

;  close,3

endfor

closedevice

save, time_save, cpcp_save, north_save, f13_save, file = 'pot.save'

end
