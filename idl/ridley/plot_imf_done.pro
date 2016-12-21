
common ffinfo, header

dt = 1.5*60.0*60.0

min_slope  = 1.0
min_slopez = 1.0
missing = -1.0e32
re = 6371.2

infile = ''
read, 'Enter input file name : ',infile

openr,1, infile

pos = fltarr(5,4)

ych = float(!d.y_ch_size)/float(!d.y_size)

ppp = 3
space = 0.01
pos_space,ppp,space,sizes, nx = 1
xmargin = 0.31
ymargin = 0.01 + 1.5*ych
for i=0,2 do begin
  get_position, ppp, space, sizes, i, tpos, 				$
                xmargin = xmargin, ymargin = ymargin, /rect
  pos(i,*) = tpos(*)
  pos(i,[0,2]) = pos(i,[0,2]) - xmargin
  pos(i,[1,3]) = pos(i,[1,3]) - ymargin
endfor

ppp = 2
space = ych*2.0*2.0
pos_space,ppp,space,sizes, nx = 1
xmargin = 0.75
ymargin = 0.00
for i=3,4 do begin
  get_position, ppp, space, sizes, i-3, tpos, 				$
                xmargin = xmargin, ymargin = ymargin, /rect
  pos(i,*) = tpos(*)
endfor

itime = intarr(6)

dx=0.0
dy=0.0
dz=0.0

bx=0.0
by=0.0
bz=0.0

dt_prop = 0.0
dt_x    = 0.0
dt_p    = 0.0

dt_b    = fltarr(6)

while not eof(1) do begin

  readf,1, itime
  readf,1, dx,dy,dz,bx,by,bz
  readf,1, dt_prop,dt_x,dt_p
  readf,1, dt_b

  if (bx ne 0.0) then begin
    yslope = by/bx
    zslope = bz/bx
  endif else begin
    yslope = by/0.01
    zslope = bz/0.01
  endelse

  if (abs(yslope) lt min_slope) then begin
    if yslope eq 0.0 then yslope = 10000.0				$
    else yslope = min_slope*yslope/abs(yslope)
  endif

  if (abs(zslope) lt min_slopez) then begin
    if zslope eq 0.0 then zslope = 10000.0				$
    else zslope = min_slope*zslope/abs(zslope)
  endif

  sy = tostr(itime(0))
  sm = chopr('0'+tostr(itime(1)),2)
  sd = chopr('0'+tostr(itime(2)),2)
  sh = chopr('0'+tostr(itime(3)),2)
  si = chopr('0'+tostr(itime(4)),2)

  psfile = 'imf'+sy+sm+sd+sh+si+'.ps'
  setdevice, psfile,'l',nothing,0.9

  marked = itime(5)
  itime(5) = 0

  c_a_to_r, itime, midtime

  year = itime(0)

  stime = midtime - dt
  etime = midtime + dt

  print, 'Reading data...'

  filename = 'imf_'+chopr('0'+tostr(year),2)+'_cross'

  col_scal = [0,1,2,13,14,15,29,30,31,45,26,27,28,4,5,6,17,18,19,7]

  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,	$
	filename = filename

  nrows = nrows(0)

  data = fltarr(nrows)

  loc = where((data_scal eq missing) or (abs(data_scal) gt 1000.0), count)
  if count gt 0 then data_scal(loc) = 1000.0

  time = time - stime

  time_axis, stime, etime, srtime, ertime,        			$
        xtickname1, xtitle1, xtickvalue, xminor, xtickn
  xtickname2 = strarr(10)+' '
  xtitle2 = ' '

  mtime = midtime-stime

  plotdumb

  mx_wind = [mean(data_scal(13,*))]
  my_wind = [mean(data_scal(14,*))]
  mz_wind = [mean(data_scal(15,*))]

  mx_imp8 = [mean(data_scal(16,*))]
  my_imp8 = [mean(data_scal(17,*))]
  mz_imp8 = [mean(data_scal(18,*))]

  range = 0.0
  for col = 0, 2 do 							$
    if max(data_scal(col,*))-min(data_scal(col,*)) gt range then	$
      range = max(data_scal(col,*))-min(data_scal(col,*))

  range = 1.2*range/2.0

  for col = 0, 2 do begin

    if col lt 2 then begin
      xtitle = xtitle2
      xtickname = xtickname2
    endif else begin
      xtitle = xtitle1
      xtickname = xtickname1
    endelse

    if col eq 0 then ytitle = 'Bx (nT)'
    if col eq 1 then ytitle = 'By (nT)'
    if col eq 2 then ytitle = 'Bz (nT)'

    data(*) = data_scal(col,*)
    m = mean(mm(data))
    yrange = [m - range, m + range]

    plot, time(col,*), data,						$
	  ystyle = 1, xstyle = 1,					$
	  /noerase, pos = pos(col,*), 					$
          yrange = yrange, max_value = 20.0,				$
          xtickname = xtickname,					$
          xtickv = xtickvalue, xticks = xtickn,				$
          xminor = xminor, xtitle=xtitle,				$
	  xrange = [srtime, ertime],					$
	  ytitle = ytitle

    data(*) = data_scal(col+3,*)
    oplot, time(col+3,*), data, max_value = 20.0, linestyle = 1

    oplot, [srtime,ertime], [0.0,0.0], linestyle = 2

    top = m+range
    bot = m-range
    dy  = 2.0*range/15.0

    oplot, [mtime,mtime], [-100,100]
    xyouts, mtime+60.0, m+range-dy, 'WIND', charsize = 0.8

    oplot, [mtime+dt_prop*60.0,mtime+dt_prop*60.0], [-100,100], linestyle = 2
    xyouts, mtime+dt_prop*60.0+60.0, m-range+dy, 'IMP 8', charsize = 0.8

    oplot, [mtime+dt_x*60.0,mtime+dt_x*60.0], [-100,100], linestyle = 1
    xyouts, mtime+dt_x*60.0+60.0, bot+2*dy, 'X', charsize = 0.8

    oplot, [mtime+dt_p*60.0,mtime+dt_p*60.0], [-100,100], linestyle = 1
    xyouts, mtime+dt_p*60.0+60.0, bot+3*dy, 'P', charsize = 0.8

    oplot, [mtime+dt_b(0)*60.0,mtime+dt_b(0)*60.0], [-100,100], linestyle = 1
    xyouts, mtime+dt_b(0)*60.0+60.0, top-2*dy, 'B', charsize = 0.8

    oplot, [mtime+dt_b(1)*60.0,mtime+dt_b(1)*60.0], [-100,100], linestyle = 1
    xyouts, mtime+dt_b(1)*60.0+60.0, top-3*dy, 'T', charsize = 0.8

    if marked eq col then 						$
      xyouts, ertime-120.0, top-dy, 'Selected Component', 		$
      alignment=1.0, charsize = 0.8

  endfor

  range = 100.0

  pn = 3

  plot, [range,-range],[0,250], /nodata, pos = pos(pn,*), 		$
	/noerase, xtitle = 'Z GSM', ytitle = 'X GSM',			$
	xrange = [range,-range], xstyle = 1

; fake magnetopause x - z

  xp = [0.0,15.0,15.0,0.0]
  zp = [15.0,15.0,-15.0,-15.0]
  oplot, zp, xp

; earth x - z

  r = 1.0
  t = findgen(19)*!pi/18.0
  xp = r*cos(t)
  yp = r*sin(t)
  oplot, xp, yp

  for j=-range,range,25 do 						$
    oplot, [j,j], [0,250], linestyle = 1

  for k=0,250,25 do							$
    oplot, [-range,range], [k,k], linestyle = 1

  oplot, mz_wind, mx_wind, psym = 1, min_value = -1000.0
  oplot, mz_imp8, mx_imp8, psym = 2, min_value = -1000.0

  xyouts, mz_wind-3, mx_wind, 'Wind'
  xyouts, mz_imp8-3, mx_imp8, 'Imp 8'

  xs = [mz_wind, mz_wind + 200.0*zslope]
  ys = [mx_wind, mx_wind + 200.0]
  oplot, xs,ys, linestyle = 3
  xs = [mz_wind, mz_wind - 200.0*zslope]
  ys = [mx_wind, mx_wind - 200.0]
  oplot, xs,ys, linestyle = 3

;  xs = [my_wind, my_wind + 200.0]
;  ys = [mx_wind, mx_wind]
;  oplot, xs,ys, linestyle = 2
;  xs = [my_wind, my_wind - 200.0]
;  ys = [mx_wind, mx_wind]
;  oplot, xs,ys, linestyle = 2

  xs = [my_wind, my_wind + 200.0]
  ys = [mx_wind, mx_wind]
  oplot, xs,ys, linestyle = 2
  xs = [my_wind, my_wind - 200.0]
  ys = [mx_wind, mx_wind]
  oplot, xs,ys, linestyle = 2

  pn = 4

  plot, [range,-range],[0,250], /nodata, pos = pos(pn,*), 		$
	/noerase, xtitle = 'Y GSM', ytitle = 'X GSM',		$
	xrange = [range,-range], xstyle = 1

; earth x - y

  oplot, xp, yp

; fake magnetopause x - y

  r = 15.0
  t = findgen(19)*!pi/18.0
  xp = r*cos(t)
  yp = r*sin(t)

  oplot, xp, yp

  for j=-range,range,25 do 						$
    oplot, [j,j], [0,250], linestyle = 1

  for k=0,250,25 do							$
    oplot, [-range,range], [k,k], linestyle = 1

  oplot, my_wind, mx_wind, psym = 1, min_value = -1000.0
  oplot, my_imp8, mx_imp8, psym = 2, min_value = -1000.0

  xyouts, my_wind-3, mx_wind, 'Wind'
  xyouts, my_imp8-3, mx_imp8, 'Imp 8'

  xs = [my_wind, my_wind + 200.0*yslope]
  ys = [mx_wind, mx_wind + 200.0]
  oplot, xs,ys, linestyle = 3
  xs = [my_wind, my_wind - 200.0*yslope]
  ys = [mx_wind, mx_wind - 200.0]
  oplot, xs,ys, linestyle = 3

  xs = [my_wind, my_wind + 200.0*(-1.0)]
  ys = [mx_wind, mx_wind + 200.0]
  oplot, xs,ys, linestyle = 2
  xs = [my_wind, my_wind - 200.0*(-1.0)]
  ys = [mx_wind, mx_wind - 200.0]
  oplot, xs,ys, linestyle = 2

  device, /close

endwhile

close,/all
set_plot, 'X'

end
