pro min_slope, data_x, data_y, slope, y_int, 			$
               init_angle=init_angle,				$
               init_y_int=init_y_int,				$
               min_error = min_error,				$
               niters    = niters

  if n_elements(init_y_int) eq 0 then init_y_int = 0.0
  if n_elements(init_angle) eq 0 then init_angle = 0.0
  if n_elements(min_error) eq 0 then min_error = 1.0e32
  if n_elements(niters) eq 0 then niters = 0

  init_slope = tan(init_angle*!pi/180.0)

  new_y = data_y - init_slope*data_x - init_y_int

  error = stdev(new_y)

  if (niters lt 128) and (error lt min_error) then begin

    if niters eq 0 then begin
      new_angle_pos = 89.99
      new_angle_neg = -89.99
    endif else begin
      new_angle_pos = init_angle + 90.0/float(niters+1)
      new_angle_neg = init_angle - 90.0/float(niters+1)
    endelse

    new_slope_pos = tan(new_angle_pos*!pi/180.0)
    new_slope_neg = tan(new_angle_neg*!pi/180.0)

    y_pos = data_y - new_slope_pos*data_x - init_y_int
    y_neg = data_y - new_slope_neg*data_x - init_y_int

    error_pos = stdev(y_pos)
    error_neg = stdev(y_neg)

    if error_pos lt error_neg then begin
      init_angle = new_angle_pos
      y_int = mean(y_pos + init_y_int)
      slope = new_slope_pos
    endif else begin
      init_angle = new_angle_neg
      y_int = mean(y_neg + init_y_int)
      slope = new_slope_neg
    endelse

    init_y_int = 0.0

    niters = niters + 1

    min_slope, data_x, data_y, slope, y_int, 			$
               init_angle=init_angle,				$
               init_y_int=init_y_int,				$
               niters    = niters
;               min_error = min_error,				$

  endif

  return

end

pro max_slope, data_x, data_y, slope,	 			$
               init_angle = init_angle,				$
               niters     = niters,				$
               percent    = percent

  if n_elements(init_angle) eq 0 then init_angle = 0.0
  if n_elements(niters) eq 0 then niters = 0
  if n_elements(percent) eq 0 then percent = 0.9

  if percent gt 1.0 then percent = percent/100.0

  init_slope = tan(init_angle*!pi/180.0)

  new_y = data_y - init_slope*data_x

  loc = where(new_y lt 0.0, count)
  per = float(count)/float(n_elements(new_y))

  if (niters lt 64) then begin

    if niters eq 0 then begin
      new_angle_pos = 89.99
      new_angle_neg = -89.99
    endif else begin
      new_angle_pos = init_angle + 90.0/float(niters+1)
      new_angle_neg = init_angle - 90.0/float(niters+1)
    endelse

    new_slope_pos = tan(new_angle_pos*!pi/180.0)
    new_slope_neg = tan(new_angle_neg*!pi/180.0)

    y_pos = new_slope_pos*data_x - data_y
    y_neg = new_slope_neg*data_x - data_y

    loc = where(y_pos gt 0.0, count)
    per_pos = float(count)/float(n_elements(new_y))

    loc = where(y_neg gt 0.0, count)
    per_neg = float(count)/float(n_elements(new_y))

    if abs(per_pos-percent) lt abs(per_neg-percent) then begin
      init_angle = new_angle_pos
      slope = new_slope_pos
    endif else begin
      init_angle = new_angle_neg
      slope = new_slope_neg
    endelse

    niters = niters + 1

    max_slope, data_x, data_y, slope,	 			$
               init_angle=init_angle,				$
               niters    = niters,				$
	       percent   = percent

  endif

  return

end

pro do_slope, x, y, pos

  ych = 1.25*float(!d.y_ch_size)/float(!d.y_size)

  min_slope, x, y, slope, y_int

  max_slope, x, y, mslope, percent = 0.95

  print, slope, y_int, mslope
  oplot, [0.0,1000.0],[y_int,1000.0*slope+y_int]
  oplot, [0.0,1000.0],[0.0,1000.0*mslope], linestyle = 2

  sslope = tostrf(slope)
  sslope = strmid(sslope,0,strpos(sslope,'.')+4)+' min/R!le!X'

  syint = tostrf(y_int)
  syint = strmid(syint,0,strpos(syint,'.')+3)+' min'

  dt_10 = 10.0*slope+y_int
  dt_50 = 50.0*slope+y_int
  sdt_10 = tostrf(dt_10)
  sdt_10 = strmid(sdt_10,0,strpos(sdt_10,'.')+2)+' min'
  sdt_50 = tostrf(dt_50)
  sdt_50 = strmid(sdt_50,0,strpos(sdt_50,'.')+2)+' min'

  xyouts, pos(2)+0.05, pos(3)-1.0*ych, 'Slope = '+sslope, /norm
  xyouts, pos(2)+0.05, pos(3)-2.0*ych, 'Y intercept = '+syint, /norm
  xyouts, pos(2)+0.05, pos(3)-3.0*ych, '!9D!Xt(10) = '+sdt_10, /norm
  xyouts, pos(2)+0.05, pos(3)-4.0*ych, '!9D!Xt(50) = '+sdt_50, /norm

  sslope = tostrf(mslope)
  sslope = strmid(sslope,0,strpos(sslope,'.')+4)+' min/R!le!X'

  dt_10 = 10.0*mslope
  dt_50 = 50.0*mslope
  sdt_10 = tostrf(dt_10)
  sdt_10 = strmid(sdt_10,0,strpos(sdt_10,'.')+2)+' min'
  sdt_50 = tostrf(dt_50)
  sdt_50 = strmid(sdt_50,0,strpos(sdt_50,'.')+2)+' min'

  xyouts, pos(2)+0.05, pos(3)-6.0*ych, 'Max Slope = '+sslope, /norm
  xyouts, pos(2)+0.05, pos(3)-7.0*ych, '!9D!Xt(10) = '+sdt_10, /norm
  xyouts, pos(2)+0.05, pos(3)-8.0*ych, '!9D!Xt(50) = '+sdt_50, /norm

  y_neg = (slope*x+y_int) - y
  loc   = where(y_neg gt 0.0, count)
  per   = 100.0*float(count)/float(n_elements(y)) + 0.5

  y_neg = (mslope*x) - y
  loc   = where(y_neg gt 0.0, count)
  mper  = 100.0*float(count)/float(n_elements(y)) + 0.5

  xyouts, pos(2)+0.05, pos(3)-10.0*ych, 			$
          '% below Slope = '+tostr(per), /norm
  xyouts, pos(2)+0.05, pos(3)-11.0*ych, 			$
          '% below Max Slope = '+tostr(mper), /norm

  return

end


;---------------------------------------------------------------------------
; start of main program

openr,1,'times_corrected.dat'

itime = intarr(6)

dx=0.0
dy=0.0
dz=0.0

dist = fltarr(1000,4)

bx=0.0
by=0.0
bz=0.0

field = fltarr(1000,4)

dt_prop = 0.0
dt_x    = 0.0
dt_p    = 0.0

dt_b    = fltarr(6)

times   = fltarr(1000,9)

n = 0

while not eof(1) do begin

  readf,1, itime
  readf,1, dx,dy,dz,bx,by,bz
  readf,1, dt_prop,dt_x,dt_p
  readf,1, dt_b

  dist(n,0:2)  = [dx,dy,dz]
  dist(n,3)    = (dx^2.0+dy^2.0+dz^2.0)^0.5

  field(n,0:2) = [bx,by,bz]
  field(n,3)   = (bx^2.0+by^2.0+bz^2.0)^0.5

  times(n,0:2) = [dt_prop,dt_x,dt_p]
  times(n,3:8) = dt_b

  n = n + 1

endwhile

close,1

dist  = dist(0:n-1,0:3)
field = field(0:n-1,0:3)
times = times(0:n-1,0:8)
times_dt = times

for i=0,8 do times_dt(*,i) = times_dt(*,i) - times(*,0)

dxy = (dist(*,2)^2.0+dist(*,1)^2.0)^0.5

setdevice, 'imf1.ps','p'

!p.charsize=1.25

ppp = 1
space = 0.15
pos_space,ppp,space,sizes
get_position, ppp, space, sizes, 0, pos
pos([0,2]) = pos([0,2]) - 0.1

yrange = mm(abs(times_dt))
xrange = mm(abs(dxy))

plot, dxy, abs(times_dt(*,1)), pos = pos, psym = 4,			$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = '!9D!XT (X distance)',					$
      yrange = yrange, xrange = xrange

do_slope, dxy, abs(times_dt(*,1)), pos

device, /close
setdevice, 'imf2.ps','p'

prompt_for_next
plot, dxy, abs(times_dt(*,2)), pos = pos, psym = 4, 			$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = '!9D!XT (Parker spiral)',				$
      yrange = yrange, xrange = xrange

do_slope, dxy, abs(times_dt(*,2)), pos

device, /close
setdevice, 'imf3.ps','p'

prompt_for_next
plot, dxy, abs(times_dt(*,3)), pos = pos, psym = 4, 			$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = '!9D!XT (Bxy Field)',					$
      yrange = yrange, xrange = xrange

do_slope, dxy, abs(times_dt(*,3)), pos

device, /close
setdevice, 'imf4.ps','p'

prompt_for_next
plot, dxy, abs(times_dt(*,4)), pos = pos, psym = 4,			$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = 'Delta T (Bt Field)',					$
      yrange = yrange, xrange = xrange

do_slope, dxy, abs(times_dt(*,4)), pos

device, /close
setdevice, 'imf5.ps','p'

prompt_for_next
loc = where(abs(dist(*,1)) le 5.0)
plot, dxy(loc), abs(times_dt(loc,1)), pos = pos, psym = 4,		$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = 'Delta T (X distance)',					$
      yrange = yrange, xrange = xrange,					$
      title = 'abs(!7D!XY) < 5.0'

do_slope, dxy(loc), abs(times_dt(loc,1)), pos

device, /close
setdevice, 'imf6.ps','p'

prompt_for_next
loc = where(abs(dist(*,2)) le 5.0)
plot, dxy(loc), abs(times_dt(loc,1)), pos = pos, psym = 4,		$
      xtitle = 'Distance in Z-Y plane',					$
      ytitle = 'Delta T (X distance)',					$
      yrange = yrange, xrange = xrange,					$
      title = 'abs(!9D!XZ) < 5.0'

do_slope, dxy(loc), abs(times_dt(loc,1)), pos

device, /close
setdevice, 'imf7.ps','p'

prompt_for_next
loc = where(dxy le 10.0)
plot, abs(dist(loc,0)), abs(times_dt(loc,1)), pos = pos, psym = 4,	$
      xtitle = 'Distance in X direction',				$
      ytitle = 'Delta T (X distance)',					$
      yrange = yrange,							$
      title = '!9D!XYZ < 10.0'

do_slope, abs(dist(loc,0)), abs(times_dt(loc,1)), pos

device, /close
set_plot, 'X'

end