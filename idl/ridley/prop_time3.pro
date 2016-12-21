pro prop_time3, pos, b, v, rho, pos2, title = title, psfile = psfile

  miu    = 4.0*!pi*1.0e-7           ; whatever
  mp     = 1.62726e-27              ; kg
  bearth = 30319.0e-9               ; Tesla
  r0     = 6371.2                   ; km
  fac    = 2.0*0.85

  if n_elements(pos) eq 0 then begin

    print, 'In order to run this program, you need to input values.'
    print, 'For a minimum run, you must input the following...'
    print, ' '

    posx = 0.0
    posy = 0.0
    print, 'Location of the satellite (1st if more than one) (x,y) :'
    read, posx, posy
    pos = [posx,posy,0]

    bx = 0.0
    by = 0.0
    print, 'IMF Bx and By values (if you don''t have these, enter 0.0,0.0):'
    read, bx,by
    b = [bx,by,0]

    vx = 0.0
    print, 'Solar wind speed (Guess if you don''t have information ~400) :'
    read, vx
    v = [-vx,0,0]

    rho = 0.0
    print, 'Solar wind density (Guess if you don''t have information ~5) :'
    read, rho

    p2x = 0.0
    p2y = 0.0
    print, 'You can also have a second satellite. If you would like to put'
    print, 'a second satellite in, enter it''s position. If not, simply
    print, 'enter 0.0,0.0 :'
    read, p2x,p2y
    if (p2x^2+p2y^2) gt 0 then pos2 = [p2x,p2y,0.0]

    title = ' '
    print, 'You can also enter a title (press return for no title) :'
    read, title

    psfile = ''
    print, 'Enter a PS file to put plot :'
    read, psfile

  endif

  if n_elements(psfile) gt 0 then setdevice, psfile, 'l', 4, 0.95

  if n_elements(v) eq 0 then v = [-400.0,0.0,0.0]
  if n_elements(rho) eq 0 then rho = 5
  if n_elements(title) eq 0 then title = ''

  x = pos(0)
  y = pos(1)
  z = pos(2)

  if sqrt(x^2.0+y^2.0+y^2.0) lt r0 then begin
    x = x*r0
    y = y*r0
    z = z*r0
  endif

  bx = b(0)
  by = b(1)
  bz = b(2)

  vx = v(0)
  vy = v(1)
  vz = v(2)

  if n_elements(pos2) gt 0 then begin
    x2 = pos2(0)
    y2 = pos2(1)
    z2 = pos2(2)
    s2 = 1
    if sqrt(x2^2.0+y2^2.0+z2^2.0) lt r0 then begin
      x2 = x2*r0
      y2 = y2*r0
      z2 = z2*r0
    endif
  endif else s2 = 0

  vmag = sqrt(vx^2.0+vy^2.0+vz^2.0)

  x_mp = r0*((2.0*bearth)^2.0/				$
	 (fac*miu*mp*(rho*1.0e6)*(vmag*1.0e3)^2.0))^(1.0/6.0)

; x_bs is the solution to the following equation:
;  x = 0.5*at^2
;  a = deceleration of the solar wind through the bowshock
;    = (v_magnetosheath - 0)/t
;    ... v_magnetosheath = 1/3 v_solarwind (spreiter and stahara, 1966)
;  t = time through the bowshock
;  so....
;  x = (1/2)*(1/3 v_sw)*t
;  6x/v_sw = t
;    we want a dt, since we are going to add this time to the time to the
;      magnetopause (with no bowshock), so (taking away 1x)
;  dt = 5x/v_sw
;  so, we can say that x_bs = 5.0*x_bowshock, so in our final equation
;      we will have somthing like :
;  t = (x_mp + x_bs)/v_sw

;    if we assume that x_bowshock = x_mp/3.0

  bowshock = x_mp/3.0
  x_bs = 5.0*bowshock

  xmax = max([x*1.2/r0,30])
  xmin = -10.0
  ymax = max([abs(y)*1.2/r0,30])
  ymin = -ymax

  xmax = 50.0
  xmin = -20.0
  ymax = 35.0
  ymin = -35.0

  pos_space,1,0.15,sizes
  get_position,1,0.05,sizes,0,plotpos
  dx = plotpos(2) - plotpos(0)
  dy = plotpos(3) - plotpos(1)
  plotpos = [0.05,0.05,0.05+dx,0.05+dy]
  xposm = 0.05+dx+0.02

  plot, [xmin,xmax],[ymin,ymax], /nodata, xstyle=1, ystyle=1, 	$
	pos = plotpos, xtitle = 'Re', 		$
	ytitle = 'Re', title = title, charsize=1.2

  a = findgen(16)*2.0*!pi/15
  usersym, cos(a), sin(a)

  oplot, [0.0], [0.0], psym = 8
  xyouts, 1.0, 1.0, 'Earth', charsize=1.2

  oplot, [x/r0], [y/r0], psym = 4
  xyouts, x/r0+1.0, y/r0+1.0, 'WIND'

  if s2 then begin
    oplot, [x2/r0], [y2/r0], psym = 4
    xyouts, x2/r0+1.0, y2/r0+1.0, 'IMP 8'
    x3 = 0.0
    y3 = -6.6*r0
    oplot, [x3/r0], [y3/r0], psym = 4
    xyouts, x3/r0+1.0, y3/r0+1.0, 'GOES 7'
  endif

  for i = xmin,xmax,10 do oplot, [i,i],[ymin,ymax],linestyle=1
  for j = 0,ymax,10 do oplot, [xmin,xmax],[j,j],linestyle=1
  for j = -10,ymin,-10 do oplot, [xmin,xmax],[j,j],linestyle=1

  xl = findgen(300)/8.0
  xl = xl - max(xl) + x_mp/r0
  yl = 4.0*sqrt(x_mp/r0-xl)
  oplot, xl,yl
  oplot, xl,-yl
  xl = xl + 2.0
  yl = 6.2*sqrt(x_mp/r0+2.0-xl)
  oplot, xl,yl
  oplot, xl,-yl

  xr = xmax - xmin
  yr = ymax - ymin

  if y gt 0 then yout = ymin + yr/20.0 else yout = ymax - yr/10.0
  oplot, [xmin/2 + xr/20.0, xmin/2 + xr/20.0 + xr/15.0],		$
	 [yout, yout]

  xyouts, xmin/2 + xr/20.0 + xr/14.0, yout, 'Propagation Plane', charsize=1.0

;  oplot, [xmin/2 + xr/4.0, xmin/2 + xr/4.0 + xr/15.0],		$
;	 [yout, yout], linestyle=2

;  xyouts, xmin/2 + xr/4.0 + xr/14.0, yout, 'Parker Spiral', charsize=1.2

  if bx ne 0 then slope = by/bx else slope = by/0.001
  if slope eq 0 then slope = 0.001
  if (bx^2.0 + by^2.0) eq 0.0 then slope = -1.0

  ys = slope*(xmin) - slope*(x/r0) + y/r0
  ye = slope*(xmax) - slope*(x/r0) + y/r0
  oplot, [xmin,xmax],[ys,ye]

  psslope = -1.0
;  ys = psslope*(xmin) - psslope*(x/r0) + y/r0
;  ye = psslope*(xmax) - psslope*(x/r0) + y/r0
;  oplot, [xmin,xmax],[ys,ye], linestyle = 2

  x_b = r0*(0.0 - y/r0 + slope*(x/r0))/slope
  psx_b = r0*(0.0 - y/r0 + psslope*(x/r0))/psslope

  x_only_time = (x-x_mp)/(-1.0*vx)
  x_mag_time = (x_b-x_mp)/(-1.0*vx)
  x_mag_bs_time = (x_b-x_mp+x_bs)/(-1.0*vx)
  x_bs_time = (x-x_mp+x_bs)/(-1.0*vx)

  psx_mag_time = (psx_b-x_mp)/(-1.0*vx)
  psx_mag_bs_time = (psx_b-x_mp+x_bs)/(-1.0*vx)

  if s2 then begin
    s2x_b = r0*(y2/r0 - y/r0 + slope*(x/r0))/slope
    s2psx_b = r0*(y2/r0 - y/r0 + psslope*(x/r0))/psslope
    s2x_only_time = (x-x2)/(-1.0*vx)
    s2x_mag_time = (s2x_b-x2)/(-1.0*vx)
    s2psx_mag_time = (s2psx_b-x2)/(-1.0*vx)
  endif

  xyouts, xposm, 0.05+dy, 'Travel Times',/norm, charsize=1.2

;  xyouts, 0.76, 0.85, 'X distance only',/norm, charsize=1.2
;  mm = fix(x_only_time/60.0)
;  ss = abs(fix(x_only_time-mm*60))
;  smm = chopn(mm,2)+' min '
;  sss = chopr('0'+tostr(ss),2)
;  if mm eq 0 then smm = ''
;  if mm eq 0 and x_only_time lt 0 then smm = '-'
;  xyouts, 0.78, 0.82, 'Subsolar : '+smm+sss+' sec',/norm
  if s2 then begin
    mm = fix(s2x_only_time/60.0)
    ss = abs(fix(s2x_only_time-mm*60))
    smm = chopn(mm,2)+' min '
    sss = chopr('0'+tostr(ss),2)
    if mm eq 0 then smm = ''
    if mm eq 0 and s2x_only_time lt 0 then smm = '-'
;    xyouts, 0.78, 0.79, 'S/C #2 : '+smm+sss+' sec',/norm
  endif

;  xyouts, 0.76, 0.75, 'X distance + Bow Shock',/norm, charsize=1.2
  mm = fix(x_bs_time/60.0)
  ss = abs(fix(x_bs_time-mm*60))
  smm = chopn(mm,2)+' min '
  sss = chopr('0'+tostr(ss),2)
  if mm eq 0 then smm = ''
  if mm eq 0 and x_bs_time lt 0 then smm = '-'
;  xyouts, 0.78, 0.72, 'Subsolar : '+smm+sss+' sec',/norm

;  xyouts, 0.76, 0.65, 'Mag field to Y location',/norm, charsize=1.2
  mm = fix(x_mag_time/60.0)
  ss = abs(fix(x_mag_time-mm*60))
  smm = chopn(mm,2)+' min '
  sss = chopr('0'+tostr(ss),2)
  if mm eq 0 then smm = ''
  if mm eq 0 and x_mag_time lt 0 then smm = '-'
;  xyouts, 0.78, 0.62, 'Subsolar : '+smm+sss+' sec',/norm
  if s2 then begin
    mm = fix(s2x_mag_time/60.0)
    ss = abs(fix(s2x_mag_time-mm*60))
    smm = chopn(mm,2)+' min '
    sss = chopr('0'+tostr(ss),2)
    if mm eq 0 then smm = ''
    if mm eq 0 and s2x_mag_time lt 0 then smm = '-'
    xyouts, xposm, 0.05+dy-0.05, 'WIND to IMP 8 : '+smm+sss+' sec',/norm
  endif

;  xyouts, 0.76, 0.55, 'Mag field + Bow Shock',/norm, charsize=1.2
  mm = fix(x_mag_bs_time/60.0)
  ss = abs(fix(x_mag_bs_time-mm*60))
  smm = chopn(mm,2)+' min '
  sss = chopr('0'+tostr(ss),2)
  if mm eq 0 then smm = ''
  if mm eq 0 and x_mag_bs_time lt 0 then smm = '-'
  xyouts, xposm, 0.05+dy-0.1, 'WIND to Subsolar : '+smm+sss+' sec',/norm
  xyouts, xposm, 0.05+dy-0.15, '(Bowshock is included)',/norm

;  xyouts, 0.76, 0.45, 'Parker Spiral to Y location',/norm, charsize=1.2
  mm = fix(psx_mag_time/60.0)
  ss = abs(fix(psx_mag_time-mm*60))
  smm = chopn(mm,2)+' min '
  sss = chopr('0'+tostr(ss),2)
  if mm eq 0 then smm = ''
  if mm eq 0 and psx_mag_time lt 0 then smm = '-'
;  xyouts, 0.78, 0.42, 'Subsolar : '+smm+sss+' sec',/norm
  if s2 then begin
    mm = fix(s2psx_mag_time/60.0)
    ss = abs(fix(s2psx_mag_time-mm*60))
    smm = chopn(mm,2)+' min '
    sss = chopr('0'+tostr(ss),2)
    if mm eq 0 then smm = ''
    if mm eq 0 and s2psx_mag_time lt 0 then smm = '-'
;    xyouts, 0.78, 0.39, 'S/C # 2 : '+smm+sss+' sec',/norm
  endif

;  xyouts, 0.76, 0.35, 'Parker Spiral + Bow Shock',/norm, charsize=1.2
  mm = fix(psx_mag_bs_time/60.0)
  ss = abs(fix(psx_mag_bs_time-mm*60))
  smm = chopn(mm,2)+' min '
  sss = chopr('0'+tostr(ss),2)
  if mm eq 0 then smm = ''
  if mm eq 0 and psx_mag_bs_time lt 0 then smm = '-'
;  xyouts, 0.78, 0.32, 'Subsolar : '+smm+sss+' sec',/norm

  imf = strcompress('('+strmid(tostrf(bx),0,5)+','+		$
			strmid(tostrf(by),0,5)+')',/remove_all)
  xyouts, xposm, 0.25, 'PP (x, y) : '+imf, /norm, charsize=1.2

  vel = strcompress(strmid(tostrf(vx),0,6),/remove_all)+' km/s'
  xyouts, xposm, 0.20, 'Speed : '+vel, /norm, charsize=1.2

  srho = strmid(tostrf(rho),0,4)+' /cm^3'
  xyouts, xposm, 0.15, 'Rho : '+srho, /norm, charsize=1.2

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end
