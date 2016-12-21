
common ffinfo, header

dt = 2.0*60.0*60.0

min_slope = 1.0

missing = -1.0e32

re = 6371.2

infile = ''
read, 'Enter input file name : ',infile

openr,1, infile

itime = intarr(6)

psfile = ''
read, 'Enter ps file name (return for screen) : ',psfile

if strlen(psfile) gt 0 then setdevice, psfile,'l',4			$
else window,0,xsize=1000,ysize=800

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
space = ych*1.75*2.0
pos_space,ppp,space,sizes, nx = 1
xmargin = 0.75
ymargin = 0.00
for i=3,4 do begin
  get_position, ppp, space, sizes, i-3, tpos, 				$
                xmargin = xmargin, ymargin = ymargin, /rect
  pos(i,*) = tpos(*)
endfor

save_cross = fltarr(1000)
save_dx    = fltarr(1000)
save_dy    = fltarr(1000)
save_dz    = fltarr(1000)
save_dt    = fltarr(1000)
save_dtx   = fltarr(1000)
save_dtb   = fltarr(1000)
save_dtt   = fltarr(1000)
save_dtp   = fltarr(1000)
npts       = 0

while not eof(1) do begin

  readf,1,itime

  marked = itime(5)
  itime(5) = 0

  c_a_to_r, itime, midtime

  year = itime(0)
  if year gt 1900 then year = year - 1900
  if year gt 100 then year = year - 100

  stime = midtime - dt
  etime = midtime + dt

  print, 'Reading data...'

  filename = 'imf_'+chopr('0'+tostr(year),2)+'_cross'

  col_scal = [0,1,2,13,14,15,29,30,31,45,26,27,28,4,5,6,17,18,19,7]

  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,	$
	filename = filename

  col_vect = [46,47,48]

  read_flat_vector, stime, etime, col_vect, time_v, data_vect, nrows,	$
	filename = filename

  nrows = nrows(0)

  if (nrows gt 2) then begin

    data = fltarr(nrows)

    loc = where((data_scal eq missing) or (abs(data_scal) gt 1000.0), count)
    if count gt 0 then data_scal(loc) = 1000.0

    time = time - stime

    time_axis, stime, etime, srtime, ertime,        			$
        xtickname1, xtitle1, xtickvalue, xminor, xtickn
    xtickname2 = strarr(10)+' '
    xtitle2 = ' '

    mtime = midtime-stime

    loc = where(time(0,*) ge mtime-120.0*5.0 and 			$
	        time(0,*) le mtime+120.0*5.0,count)

    if count gt 0 then begin

      mini = 100000.0
      j = -1

      for i=0,count-1 do begin

        loc2 = where(data_vect(0,loc(i),*) ne missing,count2)

	if count2 gt 0 then begin

	  if (abs(time(0,loc(i))-mtime) lt mini) then begin
            j = i
	    mini = abs(time(0,loc(i))-mtime)
          endif

	endif

      endfor

      if j gt 0 then begin
        loc(0) = loc(j)
        loc2 = where(data_vect(0,loc(0),*) ne missing,count2)
      endif else count2 = 0

    endif else count2 = 0

    if (count2 gt 0) then begin

      print, itime

      plotdumb

      mx_wind = [mean(data_scal(13,*))]
      my_wind = [mean(data_scal(14,*))]
      mz_wind = [mean(data_scal(15,*))]

      mx_imp8 = [mean(data_scal(16,*))]
      my_imp8 = [mean(data_scal(17,*))]
      mz_imp8 = [mean(data_scal(18,*))]

      for col = 0, 2 do begin

        if col lt 2 then begin
          xtitle = xtitle2
          xtickname = xtickname2
        endif else begin
          xtitle = xtitle1
          xtickname = xtickname1
        endelse

        data(*) = data_scal(col,*)
        plot, time(col,*), data,					$
	  ystyle = 1, xstyle = 1,					$
	  /noerase, pos = pos(col,*), 					$
          yrange = [-20,20], max_value = 20.0,				$
          xtickname = xtickname,					$
          xtickv = xtickvalue, xticks = xtickn,				$
          xminor = xminor, xtitle=xtitle,				$
	  xrange = [srtime, ertime]

        data(*) = data_scal(col+3,*)
        oplot, time(col+3,*), data, max_value = 20.0, linestyle = 1

        oplot, [srtime,ertime], [0.0,0.0], linestyle = 2

        oplot, [mtime+data_scal(6,loc(0))*60.0,			$
	        mtime+data_scal(6,loc(0))*60.0], 			$
	       [-100,100], linestyle = 3
        xyouts, mtime+data_scal(6,loc(0))*60.0+60.0, -19, 'X'

        oplot, [mtime+data_scal(8,loc(0))*60.0,			$
	        mtime+data_scal(8,loc(0))*60.0], 			$
	       [-100,100], linestyle = 3
        xyouts, mtime+data_scal(8,loc(0))*60.0+60.0, -19, 'P'

; Use the magnetic field in the X-Y plane only:

;	dt_mag = data_scal(7,loc(0))
;        oplot, [mtime+dt_mag*60.0,				$
;	        mtime+dt_mag*60.0], 				$
;	       [-100,100], linestyle = 3
;        xyouts, mtime+dt_mag*60.0+60.0, -19, 'B'

; Use the magnetic field in the X-Y plane only (10 minutes before the event):

	bx = mean(data_scal(0,loc(0)-5:loc(0)),missing)
	by = mean(data_scal(1,loc(0)-5:loc(0)),missing)
	bz = mean(data_scal(2,loc(0)-5:loc(0)),missing)

	vx = data_scal(19,loc(0)-5)

        dx = (mx_wind - mx_imp8)*re
        dy = (my_wind - my_imp8)*re
        dz = (mz_wind - mz_imp8)*re

        if (bx ne 0.0) then slope = by/bx				$
        else slope = by/0.01
 
	if (abs(slope) lt min_slope) then begin
	  if slope eq 0.0 then slope = 10000.0				$
          else slope = min_slope*slope/abs(slope)
        endif

        if (slope ne 0.0) and (bx ne missing) then begin
          magx = (-1.0*dy + slope*mx_wind*re)/slope
          dt_mag = ((magx-mx_imp8*re)/(-1.0*vx))/60.0
          oplot, [mtime+dt_mag*60.0,				$
	        mtime+dt_mag*60.0], 				$
	       [-100,100], linestyle = 3
          xyouts, mtime+dt_mag*60.0+60.0, -19, 'B'
        endif

; Use the total magnetic field (10 minutes before the event):

        if (bx ne 0.0) then begin
          yslope = by/bx
          zslope = bz/bx
        endif else begin
          yslope = by/0.01
          zslope = bz/0.01
        endelse

        if (abs(yslope) gt min_slope) and (bx ne missing) then		$
          ymagx = (-1.0*dy + yslope*mx_wind*re)/yslope			$
	else begin
	  if yslope eq 0.0 then yslope = 10000.0			$
          else yslope = min_slope*yslope/abs(yslope)
          ymagx = (-1.0*dy + yslope*mx_wind*re)/yslope
        endelse

        if (abs(zslope) gt min_slope) and (bx ne missing) then		$
          magx = (-1.0*dz + zslope*ymagx)/zslope			$
        else begin
          magx = ymagx
          zslope = 100000.0
        endelse

        if (bx ne missing) then begin
          dt_tmag = ((magx-mx_imp8*re)/(-1.0*vx))/60.0
          oplot, [mtime+dt_tmag*60.0,					$
	        mtime+dt_tmag*60.0], 					$
	       [-100,100], linestyle = 3
          xyouts, mtime+dt_tmag*60.0+60.0, -19, 'Bt'
        endif

; cross correlation time delay:

	maxc = max(data_vect(col,loc(0),loc2))
        loc3 = where(data_vect(col,loc(0),*) eq maxc) 
        t = 2.0*(loc3(0)-30.0)*60.0

        oplot, [mtime+data_scal(9,loc(0))*60.0+t,			$
	        mtime+data_scal(9,loc(0))*60.0+t], 			$
	       [-100,100], linestyle = 2
        xyouts, mtime+data_scal(9,loc(0))*60.0+t+60.0, -19, 'C'
        xyouts, mtime+data_scal(9,loc(0))*60.0+t+60.0, 18, 		$
		tostr(fix(100*maxc))

        if (col eq marked) then begin

          oplot, [mtime,mtime], [-100,100]

          save_cross(npts) = maxc
          save_dx(npts)    = data_scal(10,loc(0))
          save_dy(npts)    = data_scal(11,loc(0))
          save_dz(npts)    = data_scal(12,loc(0))
          save_dt(npts)    = data_scal(9,loc(0))*60.0 + t
          save_dtx(npts)   = data_scal(6,loc(0))*60.0
          save_dtb(npts)   = dt_mag*60.0
          save_dtt(npts)   = dt_tmag*60.0
          save_dtp(npts)   = data_scal(8,loc(0))*60.0
          npts             = npts + 1

	  xyouts, 0.0, 1.0-ych, 'Dt(X) = '+				$
              tostr(abs(save_dt(npts-1)-save_dtx(npts-1))/60.0)+' min',	$
		  /norm, charsize=1.0

	  xyouts, 0.13, 1.0-ych, 'Dt(B) = '+				$
              tostr(abs(save_dt(npts-1)-save_dtb(npts-1))/60.0)+' min',	$
		  /norm, charsize=1.0

	  xyouts, 0.26, 1.0-ych, 'Dt(Bt) = '+				$
              tostr(abs(save_dt(npts-1)-save_dtt(npts-1))/60.0)+' min',	$
		  /norm, charsize=1.0

	  xyouts, 0.39, 1.0-ych, 'Dt(P) = '+				$
              tostr(abs(save_dt(npts-1)-save_dtp(npts-1))/60.0)+' min',	$
		  /norm, charsize=1.0

	  xyouts, 0.52, 1.0-ych, 'Cross = '+				$
              tostr(100*(save_cross(npts-1)))+' %',	$
		  /norm, charsize=1.0

        endif

      endfor

      range = max([mz_wind,my_wind,mz_imp8,my_imp8,[50]])
      range = float(fix(range/50.0))*50.0 + 50.0
      range = 100.0

      pn = 3

      plot, [range,-range],[0,250], /nodata, pos = pos(pn,*), 		$
	    /noerase, xtitle = 'Z GSM', ytitle = 'X GZM',		$
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

      xs = [my_wind, my_wind + 200.0]
      ys = [mx_wind, mx_wind]
      oplot, xs,ys, linestyle = 2
      xs = [my_wind, my_wind - 200.0]
      ys = [mx_wind, mx_wind]
      oplot, xs,ys, linestyle = 2

      xs = [my_wind, my_wind + 200.0]
      ys = [mx_wind, mx_wind]
      oplot, xs,ys, linestyle = 1
      xs = [my_wind, my_wind - 200.0]
      ys = [mx_wind, mx_wind]
      oplot, xs,ys, linestyle = 1

      pn = 4

      plot, [range,-range],[0,250], /nodata, pos = pos(pn,*), 		$
	    /noerase, xtitle = 'Y GSM', ytitle = 'X GZM',		$
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

      xs = [my_wind, my_wind + 200.0*slope]
      ys = [mx_wind, mx_wind + 200.0]
      oplot, xs,ys, linestyle = 2
      xs = [my_wind, my_wind - 200.0*slope]
      ys = [mx_wind, mx_wind - 200.0]
      oplot, xs,ys, linestyle = 2

      xs = [my_wind, my_wind + 200.0*(-1.0)]
      ys = [mx_wind, mx_wind + 200.0]
      oplot, xs,ys, linestyle = 1
      xs = [my_wind, my_wind - 200.0*(-1.0)]
      ys = [mx_wind, mx_wind - 200.0]
      oplot, xs,ys, linestyle = 1

      if !d.name eq 'X' then prompt_for_next

    endif

  endif

endwhile

close,1

save_cross = save_cross(0:npts-1)
save_dx    = save_dx(0:npts-1)
save_dy    = save_dy(0:npts-1)
save_dz    = save_dz(0:npts-1)
save_dt    = save_dt(0:npts-1)
save_dtx   = abs(save_dtx(0:npts-1) - save_dt)/60.0
save_dtb   = abs(save_dtb(0:npts-1) - save_dt)/60.0
save_dtt   = abs(save_dtt(0:npts-1) - save_dt)/60.0
save_dtp   = abs(save_dtp(0:npts-1) - save_dt)/60.0

save_yz  = (save_dy^2.0 + save_dz^2.0)^0.5
save_xyz = (save_dx^2.0 + save_dy^2.0 + save_dz^2.0)^0.5

ppp = 1
space = 0.15
pos_space,ppp,space,sizes
get_position, ppp, space, sizes, 0, pos

plot, save_cross, save_yz, psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'YZ Distance'
plot, save_cross, save_xyz, psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Total Distance'
plot, save_cross, abs(save_dx), psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'X Distance'
plot, save_cross, abs(save_dy), psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Y Distance'
plot, save_cross, abs(save_dz), psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Z Distance'

plot, save_cross, save_dtx, psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Delta T (X-Distance)'
plot, save_cross, save_dtb, psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Delta T (Mag. Field)'
plot, save_cross, save_dtp, psym = 1, pos = pos, 			$
      xtitle = 'Cross Correlation', ytitle = 'Delta T (Parker)'

for cross = 100, 99 do begin

  min_cross = float(cross)/100.0

  loc = where(save_cross ge min_cross,count)

  if count gt 0 then begin

    save_cross = save_cross(loc)
    save_dx    = save_dx(loc)
    save_dy    = save_dy(loc)
    save_dz    = save_dz(loc)
    save_dt    = save_dt(loc)
    save_dtx   = save_dtx(loc)
    save_dtb   = save_dtb(loc)
    save_dtt   = save_dtt(loc)
    save_dtp   = save_dtp(loc)
    save_yz    = save_yz(loc)
    save_xyz   = save_xyz(loc)

    plot, save_yz, save_dtx, psym = 1, pos = pos, 			$
        xtitle = 'YZ Distance', ytitle = 'Delta T (X-Distance)',	$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, save_xyz, save_dtx, psym = 1, pos = pos, 			$
        xtitle = 'Total Distance', ytitle = 'Delta T (X-Distance)',	$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dx), save_dtx, psym = 1, pos = pos, 			$
        xtitle = 'X Distance', ytitle = 'Delta T (X-Distance)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dy), save_dtx, psym = 1, pos = pos, 			$
        xtitle = 'Y Distance', ytitle = 'Delta T (X-Distance)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dz), save_dtx, psym = 1, pos = pos, 			$
        xtitle = 'Z Distance', ytitle = 'Delta T (X-Distance)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)

    plot, save_yz, save_dtb, psym = 1, pos = pos, 			$
        xtitle = 'YZ Distance', ytitle = 'Delta T (Mag. Field)',	$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, save_xyz, save_dtb, psym = 1, pos = pos, 			$
        xtitle = 'Total Distance', ytitle = 'Delta T (Mag. Field)',	$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dx), save_dtb, psym = 1, pos = pos, 			$
        xtitle = 'X Distance', ytitle = 'Delta T (Mag. Field)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dy), save_dtb, psym = 1, pos = pos, 			$
        xtitle = 'Y Distance', ytitle = 'Delta T (Mag. Field)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dz), save_dtb, psym = 1, pos = pos, 			$
        xtitle = 'Z Distance', ytitle = 'Delta T (Mag. Field)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)

    plot, save_yz, save_dtp, psym = 1, pos = pos, 			$
        xtitle = 'YZ Distance', ytitle = 'Delta T (Parker)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, save_xyz, save_dtp, psym = 1, pos = pos, 			$
        xtitle = 'Total Distance', ytitle = 'Delta T (Parker)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dx), save_dtp, psym = 1, pos = pos, 			$
        xtitle = 'X Distance', ytitle = 'Delta T (Parker)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dy), save_dtp, psym = 1, pos = pos, 			$
        xtitle = 'Y Distance', ytitle = 'Delta T (Parker)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)
    plot, abs(save_dz), save_dtp, psym = 1, pos = pos, 			$
        xtitle = 'Z Distance', ytitle = 'Delta T (Parker)',		$
	title = 'Cross Correlation >= '+tostr(min_cross*100)

  endif

endfor

if !d.name eq 'PS' then begin
  device, /close
  set_plot, 'X'
endif

end
