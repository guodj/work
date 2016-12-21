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

pro do_slope, x, y, pos, mslope, slope, y_int, nolabels = nolabels

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

  if n_elements(nolabels) eq 0 then begin
    xyouts, pos(2)+0.05, pos(3)-1.0*ych, 'Slope = '+sslope, /norm
    xyouts, pos(2)+0.05, pos(3)-2.0*ych, 'Y intercept = '+syint, /norm
    xyouts, pos(2)+0.05, pos(3)-3.0*ych, '!9D!Xt(10) = '+sdt_10, /norm
    xyouts, pos(2)+0.05, pos(3)-4.0*ych, '!9D!Xt(50) = '+sdt_50, /norm
  endif

  sslope = tostrf(mslope)
  sslope = strmid(sslope,0,strpos(sslope,'.')+4)+' min/R!le!X'

  dt_10 = 10.0*mslope
  dt_50 = 50.0*mslope
  sdt_10 = tostrf(dt_10)
  sdt_10 = strmid(sdt_10,0,strpos(sdt_10,'.')+2)+' min'
  sdt_50 = tostrf(dt_50)
  sdt_50 = strmid(sdt_50,0,strpos(sdt_50,'.')+2)+' min'

  if n_elements(nolabels) eq 0 then begin
    xyouts, pos(2)+0.05, pos(3)-6.0*ych, 'Max Slope = '+sslope, /norm
    xyouts, pos(2)+0.05, pos(3)-7.0*ych, '!9D!Xt(10) = '+sdt_10, /norm
    xyouts, pos(2)+0.05, pos(3)-8.0*ych, '!9D!Xt(50) = '+sdt_50, /norm
  endif

  y_neg = (slope*x+y_int) - y
  loc   = where(y_neg gt 0.0, count)
  per   = 100.0*float(count)/float(n_elements(y)) + 0.5

  y_neg = (mslope*x) - y
  loc   = where(y_neg gt 0.0, count)
  mper  = 100.0*float(count)/float(n_elements(y)) + 0.5

  yp = 75.0
  xp = 80.0/mslope
  xyouts, xp, yp, tostr(mper)+'%', charsize = 0.8

  yp = slope*80.0 + y_int
  xp = 110.0
  
  xyouts, xp, yp, tostr(per)+'%', charsize = 0.8

  return

end


;---------------------------------------------------------------------------
; start of main program

re = 6371.2

openr,1,'times_corrected.dat'

itime = intarr(6)

dx=0.0
dy=0.0
dz=0.0

dist = fltarr(1000,4)

bx=0.0
by=0.0
bz=0.0

field = fltarr(1000,7)

dt_prop = 0.0
dt_x    = 0.0
dt_p    = 0.0

dt_b    = fltarr(6)

times   = fltarr(1000,9)

vx      = fltarr(1000)

n = 0

while not eof(1) do begin

  readf,1, itime
  readf,1, dx,dy,dz,bx,by,bz
  readf,1, dt_prop,dt_x,dt_p
  readf,1, dt_b

  print, n, itime

  dist(n,0:2)  = [dx,dy,dz]
  dist(n,3)    = (dx^2.0+dy^2.0+dz^2.0)^0.5

  field(n,0:2) = [bx,by,bz]
  field(n,3)   = (bx^2.0+by^2.0+bz^2.0)^0.5

  times(n,0:2) = [dt_prop,dt_x,dt_p]
  times(n,3:8) = dt_b

  itime(5) = 0
  c_a_to_r, itime, midtime
  year = itime(0)
  if year gt 1900 then year = year - 1900
  if year gt 100 then year = year - 100
 
  stime = midtime - 10.0*60.0
  etime = midtime + 10.0*60.0
 
  filename = 'imf_'+chopr('0'+tostr(year),2)+'_cross'
 
  col_scal = [7]
 
  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,     $
        filename = filename

  vx(n) = mean(data_scal(0,*))/re

  stime = midtime + 5.0*60.0
  etime = midtime + 15.0*60.0
 
  filename = 'imf_'+chopr('0'+tostr(year),2)+'_cross'
 
  col_scal = [1,2,3]
 
  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,     $
        filename = filename

  field(n,4) = mean(data_scal(0,*))
  field(n,5) = mean(data_scal(1,*))
  field(n,6) = mean(data_scal(2,*))

  n = n + 1

endwhile

close,1

dist  = dist(0:n-1,0:3)
field = field(0:n-1,0:6)
vx    = vx(0:n-1)
times = times(0:n-1,0:8)
times_dt = times

for i=0,8 do times_dt(*,i) = times_dt(*,i) - times(*,0)

dxy = (dist(*,2)^2.0+dist(*,1)^2.0)^0.5

a_y = 45.0
a_y_2 = 50.0
a_z = 60.0

min_slope = tan(a_y*!pi/180.0)
y_min_slope = tan(a_y_2*!pi/180.0)
z_min_slope = tan(a_z*!pi/180.0)

for i=0,n-1 do begin

  bx = field(i,0)
  by = field(i,1)
  bz = field(i,2)

  bx_a = field(i,4)
  by_a = field(i,5)
  bz_a = field(i,6)

  bb = [bx,by,bz]
  ba = [bx_a,by_a,bz_a]

  bbba = crossp(bb,ba)

  bx = bbba(1)
  by = -bbba(0)
  bz = -bbba(0)

  dx = dist(i,0)
  dy = dist(i,1)
  dz = dist(i,2)

  if (bx ne 0.0) then slope = by/bx 					$
  else slope = by/0.01
 
  if (abs(slope) lt min_slope) then begin
    if slope eq 0.0 then slope = 10000.0				$
    else slope = min_slope*slope/abs(slope)
  endif
 
  magx = (-1.0*dy + slope*dx)/slope
  dt_mag = (magx/(-1.0*vx(i)))/60.0

  times_dt(i,3) = dt_mag - times(i,0)

  bx = bbba(2)

  if (bx ne 0.0) then begin
    yslope = by/bx
    zslope = bz/bx
  endif else begin
    yslope = by/0.01
    zslope = bz/0.01
  endelse
 
  if (abs(yslope) lt y_min_slope) then begin
;      print, 'Setting yslope from ',yslope, ' to ',min_slope
    if yslope eq 0.0 then yslope = 10000.0				$
    else yslope = y_min_slope*yslope/abs(yslope)
  endif
 
  ymagx = (-1.0*dy + yslope*dx)/yslope

  if (abs(zslope) lt z_min_slope) then begin
    if zslope eq 0.0 then begin
      zslope = 100000.0
    endif else zslope = z_min_slope*zslope/abs(zslope)
  endif

  magx = (-1.0*dz + zslope*ymagx)/zslope

  dt_tmag = (magx/(-1.0*vx(i)))/60.0

  times_dt(i,4) = dt_tmag - times(i,0)

endfor

titles = ['X distance','Parker Spiral','Bxy Field','B Field',		$
	  'Bxy after','B after','Bxy cross','B cross']

;  xyouts, 0.0,0.0, 'Min slope : '+string(min_slope), charsize=0.8

yrange = mm(abs(times_dt))
xrange = mm(abs(dxy))

setdevice, 'idl.ps','p'

!p.charsize=1.25

ppp = 2
space = 0.1
pos_space,ppp,space,sizes

for j=1,8 do begin

  pn = (j-1) mod ppp
  if pn eq 0 then plotdumb

  get_position, ppp, space, sizes, pn, pos

  plot, dxy, abs(times_dt(*,j)), pos = pos, psym = 4,		$
            xtitle = 'Distance in Z-Y plane',				$
            ytitle = '!9D!XT ('+titles(j-1)+')',			$
            yrange = yrange, xrange = xrange, /noerase

  do_slope, dxy, abs(times_dt(*,j)), pos, mslope, slope, yint

endfor

closedevice

setdevice, 'imf_x.ps','p',4,0.9

!p.charsize=1.25

ppp = 4
space = 0.033
pos_space,ppp,space,sizes

loc = where(dxy lt 20.0)

for j=3,3 do begin

  if j gt 2 then begin
    xtitle =  'X Distance'
    xtickname = strarr(10)
  endif else begin
    xtitle = ''
    xtickname = strarr(10)+' '
  endelse

  if j mod 2 eq 1 then begin
    ytitle = '!9D!XT'
    ytickname = strarr(10)
  endif else begin
    ytitle = ''
    ytickname = strarr(10)+' '
  endelse

  if j eq 1 then abcd = 'A. '
  if j eq 2 then abcd = 'B. '
  if j eq 3 then abcd = ''
  if j eq 4 then abcd = 'D. '

  pn = (j-1) mod ppp
  if pn eq 0 then plotdumb

  get_position, ppp, space, sizes, pn, pos

  plot, abs(dist(loc,0)), abs(times_dt(loc,j-2)), pos = pos, psym = 4,	$
            xtitle = xtitle, xtickname = xtickname,			$
            ytitle = ytitle, ytickname = ytickname,			$
            yrange = [0,80], /noerase

  xyouts, 5.0, 72.0, abcd+titles(j-3), charsize = 1.2

;  do_slope, abs(dist(loc,0)), abs(times_dt(loc,j-2)), pos, mslope, slope, yint, /nolabels

endfor

closedevice

end
