;*****************************************************************************

pro plot_polar_vector, point, north, east, fac, pole, white = white

;*****************************************************************************

  if abs(north) lt 2000.0 and abs(east) lt 2000.0 then begin

    xn = point(0)
    yn = point(1)
    xo = pole(0)
    yo = pole(1)
;    print,'point=',xn,yn
;    print,'pole=',xo,yo		

    oplot, [xn],[yn],psym=4,color=255

    thick = 5.0

    dx = xn - xo
    dy = yn - yo
    dr = sqrt(dx^2.0+dy^2.0)
    len = sqrt(north^2.0 + east^2.0)
    if (dr gt 0.0) then begin
      theta = acos(dx/dr)
      if dy lt 0.0 then theta = -1.0*theta
      if (len gt 0) then alpha = -acos(east/len)
      if north gt 0 then alpha = -1.0*alpha
      xvec = fac*len*sin(alpha - (!pi-theta))
      yvec = -fac*len*cos(alpha - (!pi-theta))

;      xvec_a1 = 0.75*fac*len*sin(alpha - (!pi-(theta+20*!pi/180.0)))
;      yvec_a1 = -0.75*fac*len*cos(alpha - (!pi-(theta+20*!pi/180.0)))

;      xvec_a2 = 0.75*fac*len*sin(alpha - (!pi-(theta-20*!pi/180.0)))
;      yvec_a2 = -0.75*fac*len*cos(alpha - (!pi-(theta-20*!pi/180.0)))

;	  print,'north=',north
;	  print,'east=',east
;	  print,'fac=',fac
;	  print,'len=',len
;	  print,'theta=',theta    	
;	  print,'alpha=',alpha
;	  print,'dx=',dx
;	  print,'dy=',dy
;	  print,'dr=',dr
;	  print,'xvec=',xvec
;	  print,'yvec=',yvec
;          print,' '

      if n_elements(white) eq 0 then begin
        plots, [xn,xn+xvec],[yn,yn+yvec], thick=thick
;        plots, [xn+xvec,xn+xvec_a1],[yn+yvec,yn+yvec_a1], thick=thick
;        plots, [xn+xvec,xn+xvec_a2],[yn+yvec,yn+yvec_a2], thick=thick
      endif else begin
        plots, [xn,xn+xvec],[yn,yn+yvec],thick=thick,color = 255
;        plots, [xn+xvec,xn+xvec_a1],[yn+yvec,yn+yvec_a1],thick=thick,color=255
;        plots, [xn+xvec,xn+xvec_a2],[yn+yvec,yn+yvec_a2],thick=thick,color=255
      endelse
    endif

  endif

  return

end

;*****************************************************************************

pro plot_dmsp_vector, pointn, pointo, length, white=white

;*****************************************************************************

  xo = pointn(0)
  yo = pointn(1)
  xn = pointo(0)
  yn = pointo(1)

  dx = xn - xo
  dy = yn - yo
  dr = sqrt(dx^2.0+dy^2.0)
  if (dr gt 0.0) then begin
    xvec = -1.0*(-dy/dr)*length
    yvec = -1.0*(dx/dr)*length
    oplot, [xo],[yo], psym=3
    if n_elements(white) eq 0 then begin
      plots, [xo,xo+xvec],[yo,yo+yvec]
    endif else begin
      plots, [xo,xo+xvec],[yo,yo+yvec], color = 255
    endelse
  endif

  return

end

;*****************************************************************************

pro read_mag, filename, uttime, data, names

;*****************************************************************************

  openr,1,filename

  data = fltarr(4,200)
  names=sindgen(200)
;  print,'names=,',names
  flag=0
  done = 0
  n = 0

  time = 0.0
  otime=0.0 
  line = ''

  while (not done) do begin

    readf,1,line

    if eof(1) then begin
      done = 1
    endif else begin
      hh = float(strmid(line,11,2))
      mm = float(strmid(line,14,2))
      ss = float(strmid(line,17,2))
      otime=time	

	time = hh + mm/60.00 + ss/3600.00
;	print,'time=',time,'uttime=',uttime,'otime=',otime	

      if (time ge uttime) and (flag eq 0) then begin
	otime=time
	flag=1
       endif			 

      if (time eq uttime) and (flag eq 1) then begin
;        print,''
;	print,'time=',time,'uttime=',uttime,'otime=',otime	
 
	data(0,n) = float(strmid(line,28,6))
	data(1,n) = float(strmid(line,35,6))
	data(2,n) = float(strmid(line,41,12))
	data(3,n) = float(strmid(line,53,12))
	names(n)=strmid(line,24,3)


;	print,'data=',data(1,n),data(2,n),data(3,n)

        n = n + 1
;	flag=1
      endif
;	print,'n=',n
    		
      if (flag eq 1) and (time gt otime) then done = 1
    endelse

  endwhile
;	print,'names=,',names
;  print,'data=,',data(*,0:n+1) 

  data = data(*,0:n-1)
  names = names(0:n-1)
;  print,'data=,',data 
  close, 1

  return

end

;*****************************************************************************

pro read_dmsp, filename, uttime, dt, data, nskip

;*****************************************************************************

  data = fltarr(3,2.0*dt*60)

  stime = uttime*60.0 - dt
  etime = uttime*60.0 + dt

  openr,1,filename

  line = ''
  for i=1,22 do readf,1,line

  done = 0
  stime_found = 0
  n = 0

  while (not done) do begin

    readf,1,line

    if eof(1) then begin
      done = 1
    endif else begin

      hh = float(strmid(line,8,2))
      mm = float(strmid(line,10,2))
      ss = float(strmid(line,12,2))
      time = hh*60.0 + mm + ss/60.0

      if time ge stime then stime_found = 1
      if time gt etime then begin
        done = 1
      endif else begin

        if stime_found then begin

	  data(0,n) = float(strmid(line,14,9))
	  data(1,n) = float(strmid(line,23,9))
	  data(2,n) = float(strmid(line,59,9))

          n = n + 1

        endif

      endelse

      for i=1,nskip do readf,1,line

    endelse

  endwhile

  if n gt 1 then data = data(*,0:n-1) else data = data(*,0:1)

  close, 1

  return

end

;-----------------------------------------------------------------------------
; add corrected geomagnetic coordinates
;-----------------------------------------------------------------------------

;*****************************************************************************

pro add_cgm, minlat, maxlat, rotation, pos

;*****************************************************************************

  cgmgrid = '/l2/software/idl/polar/save_files/cgm_oplot_96.dat'

  openr,1,cgmgrid
  line = ''
  readf,1,line
  readf,1,line
  readf,1,line

  nlat = 90-minlat+1
  nlon = 360/5+1
  cgmdata = fltarr(4,nlat*nlon)
  bdata = fltarr(3,nlat*nlon)

  x = fltarr(18)

  for i=90,minlat,-1 do for j=0,360,5 do begin
    readf,1,line
    x(0) = strmid(line,7,6)
    x(1) = strmid(line,15,6)
    x(2) = strmid(line,21,6)
    x(3) = strmid(line,27,7)
    x(4) = strmid(line,34,7)
    x(5) = strmid(line,41,7)
    x(6) = strmid(line,48,8)
    x(7) = strmid(line,56,7)
    cgmdata(1,(90-i)*nlon+j/5) = 90.0-x(0)
    cgmdata(0,(90-i)*nlon+j/5) = float(j) + rotation
    cgmdata(2,(90-i)*nlon+j/5) = x(3)
    cgmdata(3,(90-i)*nlon+j/5) = x(4)
    bdata(*,(90-i)*nlon+j/5) = x(5:7)
  endfor

  close,1

  range = maxlat-minlat
  plot, [-range, range], [-range, range],         $
        xstyle = 5, ystyle = 5,                 $
        /noerase, pos = pos, /nodata

  levels = findgen(10)*10.0
  basic_contour_new, cgmdata(0:2,*,*), 40.0, pos, levels, /white

  return

end

;-----------------------------------------------------------------------------
; set_up_color
;-----------------------------------------------------------------------------

;*****************************************************************************

pro set_up_color, ncolors, ctname

;*****************************************************************************

  openr,11,ctname
  ncolors=0
  readf,11,ncolors
  color = fltarr(3,255)
  readf,11,color
  close,11
  tvlct,color(0,*),color(1,*),color(2,*)

  return

end

;*****************************************************************************

pro color_bar, ncolors, maxval, r

;******************************************************************************
  if r eq 0 then begin

    plot, [0,maxval], /noerase, pos = [0.97,0.3,1.0,0.7], 	$
	xstyle=5, ystyle=1, ytitle = 'photons/(cm2 s)'

  endif else begin

    plot, [-maxval,maxval], /noerase, 				$
	pos = [0.97,0.3,1.0,0.7], xstyle=5, ystyle=1,		$
	ytitle = 'photons/(cm2 s)'

  endelse

  plot, [0,ncolors], /noerase, pos = [0.97,0.3,1.0,0.7], xstyle=5,ystyle=5
  x = [0.0,0.0,1.0,1.0,0.0]
  y = [0.0,1.0,1.0,0.0,0.0]
  for i=0,ncolors-1 do 					$
    polyfill, x, float(i)+y, color = i
  plots, [0.0,0.0], [0.0,ncolors]
  plots, [1.0,1.0], [0.0,ncolors]
  plots, [0.0,1.0], [0.0,0.0]
  plots, [0.0,1.0], [ncolors,ncolors]

  return

end

;*****************************************************************************

pro polar, latitude, longitude, image, times, replot = replot

;*****************************************************************************

if n_elements(replot) eq 0 then replot = 0 else replot = 1

if replot eq 0 then begin
  done = 0
  while done eq 0 do begin
    fin = ''
    print, 'enter filename(s) to read (can use * as wildcard): '
    read, fin
    flist = findfile(fin)
    if (n_elements(flist) gt 1) then done = 1			$
    else if (strlen(flist(0)) gt 0) then done = 1			$
    else print, 'No file found, please re-enter selection'
  endwhile
  nfiles = n_elements(flist)
  print, tostr(nfiles)+' files to read'
endif else begin
  nfiles = n_elements(times)
endelse

dmspfile=''
magfile = ''

print, 'Enter file name which contains DMSP data [none default]:'
read, dmspfile
if strlen(dmspfile) gt 0 then begin
  nskip = ''
  print, 'The DMSP plots are sometimes very crowded. It may be'
  print, 'good to skip a few seconds between each vector.'
  print, 'Enter number of seconds to skip [9 default] : '
  read, nskip
  if strlen(nskip) eq 0 then nskip = 9 else nskip = fix(nskip)
endif

print, 'Enter file name which contains magnetometer data [none default]:'
read, magfile

print, '1. Plot to screen'
print, '2. Plot to PS file'
print, '3. Plot to GIF file'
print, '4. Plot only images to TIFF file'
print, 'Enter plotting location :'
plotloc = 0
read, plotloc

smaxval = ''
print, 'Enter maximum value of color table (200-750 photons/(cm2 s)) [400]:'
read, smaxval
if strlen(smaxval) gt 0 then maxval = float(smaxval) else maxval = 400.0

ctname = ''
pct = '/d/ridley/d.data/d.polar/d.software/d.idl/image.ct'
rct = '/d/ridley/d.data/d.polar/d.software/d.idl/residual.ct'
print,''
print, 'Recommended color table :'
print, '  ',pct
print, 'Recommended residual color table :'
print, '  ',rct
print, ''
print, 'Enter color table file name [default is above]:'
read, ctname

if replot ne 0 then begin
  print, 'These images may have already been smoothed. If they have, you'
  print, 'may want to answer "n" to the following question...'
endif
print, 'Would you like these images to be smoothed (y/n) ? [n]'
que = ''
read,que
if strmid(que,0,1) eq 'y' then smooth = 1 else smooth = 0

ppp = 1
rough = 0
if (plotloc le 2) then begin
  print, 'Enter number of plots per page [1 default] :'
  que = ''
  read, que
  if strlen(que) gt 0 then ppp = fix(que)
  rough = 1
  print, 'PS files to be smaller, but less quality? [y default] :'
  que = ''
  read, que
  if (strlen(que) gt 0) and (strmid(que,0,1) eq 'n') then rough = 0 
endif

if (plotloc eq 3) then begin
  print,'Enter number of pixels across that you would like.'
  print, '200 is tiny, 400 is medium, 800 is huge, default is 400.'
  que = ''
  read, que
  if strlen(que) eq 0 then gifpix=400 else gifpix=fix(que)
endif 

sub = -1
if (nfiles gt 1) then begin
  que = ''
  print, 'Would you like a subtraction to be done (y/n) ? [n]'
  read, que
  if (strmid(que,0,1) eq 'y') then begin
    print,''
    for i=0,nfiles-1 do print, tostr(i+1)+'. '+flist(i)
    print, 'Enter filename to which contains subtraction:'
    read, sub
    sub = sub - 1
    if strlen(ctname) eq 0 then ctname = rct
  endif else if strlen(ctname) eq 0 then ctname = pct
endif else if strlen(ctname) eq 0 then ctname = pct

if plotloc ne 4 then begin

  geolats = ''
  print, 'Do you want geographic latitudes and longitudes (y/n) ? [y] '
  read, geolats
  if strlen(geolats) eq 0 then geolats = 'y'

  que = ''
  print, 'Do you want corrected geomagnetic coordinates (y/n) ? [y] '
  read, que
  if strmid(que,0,1) eq 'n' then cgm = 0 else cgm = 1

  que = ''
  print, 'Do you want to zoom in on a region (y/n) ? [n]'
  read,que

  if strmid(que,0,1) eq 'y' then begin

    dlat = 25
    dlon = 100

    print, 'delta lat = ',tostr(dlat),', delta long = ',tostr(dlon)
    print, 'Select center latitude [Greenland = 70, default] : '
    que = ''
    read, que
    if strlen(que) eq 0 then clat = 70 else clat = fix(que)

    print, 'Select center longitude [Greenland = 315, default] : '
    que = ''
    read, que
    if strlen(que) eq 0 then clon = 315 else clon = fix(que)

    if clat gt 90-dlat/2 then clat = 90-dlat/2
    minlat = clat - dlat/2
    maxlat = clat + dlat/2

;    if clon gt 360-dlon/2 then clon = 360-dlon/2
;    if clon lt dlon/2 then clon = dlon/2
    minlon = clon - dlon/2
    maxlon = clon + dlon/2

  endif else begin

    minlat = 50
    maxlat = 90
    minlon = 0
    maxlon = 360

  endelse

endif else begin

  minlat = 50
  maxlat = 90
  minlon = 0
  maxlon = 360

endelse

range = maxlat-minlat
if maxlon gt 360 then begin
  minlon = minlon - 360
  maxlon = maxlon - 360
endif

limit = [minlat, minlon, maxlat, maxlon]
if ppp gt 0 then !p.multi=[0,1,ppp]

if (replot eq 0) then begin
  times     = strarr(nfiles)
  d1   = 0.0
  d2   = 0.0
  d3   = 0.0
  time = ''

  for n=0,nfiles-1 do begin

    print, 'Reading file ',flist(n)
    openr,1,flist(n)
    xsd = 0
    ysd = 0
    readf,1, xsd,ysd
    readf,1, time
    times(n) = time

    if n eq 0 then begin
      latitude  = fltarr(nfiles,xsd,ysd)
      longitude = fltarr(nfiles,xsd,ysd)
      image     = fltarr(nfiles,xsd,ysd)
      xs = xsd
      ys = ysd
    endif

    imag = fltarr(xsd,ysd)
    lat = fltarr(xsd,ysd)
    lon = fltarr(xsd,ysd)

    readu,1,lat
    readu,1,lon
    readu,1,imag

    latitude(n,*,*) = lat
    longitude(n,*,*) = lon
    image(n,*,*) = imag

    close,1

  endfor

endif

if smooth then begin
  for nn=0,nfiles-1 do begin
    print, 'Smoothing time ',times(nn)
    for i=1,xs-2 do for j=1,ys-2 do begin
      ave = 0.0
      np = 0
      for k=i-1,i+1 do for l=j-1,j+1 do begin
        if (image(nn,k,l) gt 0) then begin
          if (k+l eq 0) then begin
	    ave = ave + 4.0*image(nn,k,l)
	    np = np + 3
	  endif else ave = ave + image(nn,k,l)
	  np = np + 1
        endif
      endfor
      if (np gt 0) then image(nn,i,j) = ave/np
    endfor
  endfor
endif

; set up map region

lats = fltarr(9)
loni = fltarr(9)
i_lat = fltarr(5)
i_lon = fltarr(5)

if rough eq 1 then begin
  gifpix = 400
  pic_image = intarr(ppp,gifpix,gifpix)
  time_st = strarr(ppp)
  uttimesave = fltarr(ppp)
endif

if sub eq -1 then n1 = 0 else n1 = 1 

for nn=n1,nfiles-1 do begin

  plot_num = (nn-n1) mod ppp

  if (sub eq -1) then begin
    n = nn
    r = 0
  endif else begin
    if (nn eq 0) then begin
      n = sub
      r = 0
    endif else begin
      r = 1
      if (nn le sub) then n = nn - 1 else n = nn
    endelse
  endelse

  time = times(n)

  if (plotloc eq 1) and (nn eq 0) then window,xsize=600,ysize=600	$
  else begin
    if (rough eq 1) and (plotloc eq 2) then begin
      set_plot, 'Z'
      device, set_res=[gifpix,gifpix], Z=0
    endif
    if (plotloc eq 2) and (plot_num eq n1) and (rough eq 0) then begin
      if r eq 0 then psfile = strmid(time,0,13)+'.ps'		$
      else psfile = strmid(time,0,13)+'_sub.ps'
      setdevice, psfile, 'l',4
    endif else begin
      if (plotloc eq 3) then begin
	set_plot, 'Z'
	device, set_res=[gifpix,gifpix], Z=0
      endif
    endelse
  endelse

  space = 0.05
  if (rough eq 0) then pos_space,ppp,space,sizes			$
  else pos_space,1,0.0,sizes

  uttime = float(strmid(time,7,2)) + float(strmid(time,9,2))/60.0

  lat = 90.0
  long = 0.0
;--------------------------------------------------------------------------
  rotation = (12-uttime)*180.0/12.0 + 180.0
;  rotation = (12-uttime)*180.0/12.0 + 300

  if (rough eq 0) then get_position,ppp,space,sizes,plot_num,pos	$
  else get_position,1,space,sizes,0,pos
  pos([0,2]) = pos([0,2]) - space/2.0
  !p.region = pos
  map_set,lat, long, rotation, limit = limit, 			$
	/stereo, /noborder, /advance

  set_up_color, ncolors, ctname

  np = maxlon - minlon + 1
  la = [fltarr(np)+minlat,fltarr(np)+maxlat,minlat]
  lo = [findgen(np)+minlon,maxlon-findgen(np),minlon]
  polyfill, lo, la, color=0

  if n_elements(xs) eq 0 then begin
    xs = n_elements(longitude(0,*,0))
    ys = n_elements(longitude(0,0,*))
  endif

  for i=1,xs-2 do for j=1,ys-2 do begin

    np = 0
    tp = 0
    for k=i-1,i+1,2 do for l=j-1,j+1,2 do			$
      if (longitude(n,k,l) ge 360+minlon) then			$
	longitude(n,k,l) = longitude(n,k,l) - 360.0
    for k=i-1,i+1,2 do for l=j-1,j+1,2 do begin
      if (latitude(n,k,l) ge minlat) and			$
	 (latitude(n,k,l) le maxlat) and			$
	 (longitude(n,k,l) ge minlon) and			$
	 (longitude(n,k,l) le maxlon) then np = np + 1
      tp = tp+1
    endfor

    if tp eq np then begin

	i_lat(0) = latitude(n,i-1,j-1)
	i_lat(1) = latitude(n,i-1,j+1)
	i_lat(2) = latitude(n,i+1,j+1)
	i_lat(3) = latitude(n,i+1,j-1)
	i_lat(4) = latitude(n,i-1,j-1)
	i_lon(0) = longitude(n,i-1,j-1)
	i_lon(1) = longitude(n,i-1,j+1)
	i_lon(2) = longitude(n,i+1,j+1)
	i_lon(3) = longitude(n,i+1,j-1)
	i_lon(4) = longitude(n,i-1,j-1)
        if r eq 0 then colo = float(ncolors)*image(n,i,j)/maxval 	$
        else begin
          colo = float(ncolors)*					$
	    ((image(n,i,j)-image(sub,i,j))/2.0+(maxval/2.0))/maxval
          if colo lt 1.0 then colo = 1.0
          if colo ge 200.0 then colo = 199.0
        endelse

;---------------------------------------------------------------------------
;        polyfill, i_lon,i_lat,color=1
;---------------------------------------------------------------------------
        polyfill, i_lon, i_lat, color=colo
      endif
  endfor

  if (strmid(geolats,0,1) eq 'y') and (rough eq 0) then begin
    np = maxlon - minlon + 1
    la = [fltarr(np)+minlat,fltarr(np)+maxlat,minlat]
    lo = [findgen(np)+minlon,maxlon-findgen(np),minlon]
    for i=minlat,maxlat,10 do oplot, lo,la+i-minlat, linestyle=1, color=255
    la = [minlat,maxlat]
    lo = [0.0,0.0]
    for i=0,360,45 do oplot, lo+i,la,linestyle=1, color = 255
    MAP_CONTINENTS,MLINETHICK=1, color = 255
  endif

  ppos,cgmpos

  if strlen(dmspfile) gt 0 then begin

    dt = 15.0
    read_dmsp, dmspfile, uttime, dt, dmspdata, nskip

    for i=0,n_elements(dmspdata(0,*))-1 do begin
      if (dmspdata(0,i) ge minlat) and (dmspdata(0,i) le maxlat) and	$
	 (dmspdata(1,i) ge minlon) and (dmspdata(0,i) le maxlon) and	$
	 (dmspdata(2,i) lt 3000.0) then 	$
        dumb = convert_coord(dmspdata(1,i), dmspdata(0,i), /to_device)	$
      else begin
        dmspdata(2,i) = 0.0
        dumb = [0.0,0.0]
      endelse
      dmspdata(0,i) = dumb(0)
      dmspdata(1,i) = dumb(1)
    endfor

  endif

  if (strlen(magfile) gt 0) and (rough eq 0) then begin

    read_mag,magfile,uttime,magdata,names

; this is the geographic location of the magnetic pole for 1996 :
    glat = 81.30
    glon = 277.64	
    if minlon lt 0.0 then glon = glon - 360
    pole = convert_coord(glon, glat, /to_device)


    for i=0,n_elements(magdata(0,*))-1 do begin
    if minlon lt 0.0 then magdata(1,i)=magdata(1,i)-360.0
   
      if (magdata(0,i) ge minlat) and (magdata(0,i) le maxlat) and	$
	 (magdata(1,i) ge minlon) and (magdata(1,i) le maxlon) then 	$
        dumb = convert_coord(magdata(1,i), magdata(0,i), /to_device)	$
      else begin
        magdata(2,i) = 0.0
        dumb = [0.0,0.0]
      endelse
      magdata(0,i) = dumb(0)
      magdata(1,i) = dumb(1)
    endfor

    plot, [0,!d.x_size],[0,!d.y_size],/nodata,/noerase,		$
	xstyle=5, ystyle=5, pos = [0.0,0.0,1.0,1.0]

    north = 0.0
    east = -50.0

    loc = where((magdata(0,*) ne 0.0) and (abs(magdata(2,*)) lt 2000.0), count)
    if count gt 1 then begin

;      maxp =3.0* max(sqrt(magdata(2,loc)^2.0+magdata(3,loc)^2.0))
;      maxp = float(fix(maxp/10.0))*10.0 + 10.0
;------------------------------------------------------------------------
	maxp=200.0	
;------------------------------------------------------------------------
      fac = abs((float(!d.x_size)/5.0)/maxp)

      for i=loc(0),n_elements(magdata(0,*))-1 do begin
;      for i=loc(0)+1,n_elements(magdata(0,*))-1 do begin
;	print,'i=',i

	if magdata(2,i) ne -1.0e32 then begin
	  symbol=names(i)	
	  magpos = [magdata(0:1,i)]
	  north = magdata(3,i)
	  east = -1.0*magdata(2,i)
;	  print,'symbol=',symbol
;---------------------------------------------------------------------
;	  oplot,[magpos(0)],[magpos(1)],psym=5,color=255
;	  xyouts,magpos(0)+250.0,magpos(1)-250.0,symbol,color=255
;
;---------------------------------------------------------------------
;	print,'magpos=',magpos
;	print,'pole=',pole

          plot_polar_vector, magpos, north, east, fac, pole, /white
        endif

      endfor

    endif else begin
      fac = abs((float(!d.x_size)/5.0)/2000.0)
      maxp = 200.0
    endelse
;    print,'names=,',names
;----------------------------------------------------------------------
    pointo=[float(!d.x_size)/5.0,float(!d.y_size)/100.0+1.0]
    pointn=[float(!d.x_size)/5.0,float(!d.y_size)/100.0]
;----------------------------------------------------------------------

;    pointn=[float(!d.x_size)/5.0+1.0,float(!d.y_size)/100.0]
;    pointo=[float(!d.x_size)/5.0,float(!d.y_size)/100.0]
    length = maxp*fac/2.0
    plot_dmsp_vector, pointn, pointo, length

;    xyouts,0.21, 0.01, tostr(maxp/2.0)+' nT', /norm
;------------------------------------------------------------------------
    if rough eq 0 then begin

      xyouts,0.0,0.95, $
	strmid(time,2,2)+'-'+strmid(time,4,2)+'-'+strmid(time,0,2) ,      $
	charsize=1.2,/norm
      xyouts,0.21, 0.02,  $ 
         '      ' +tostr(maxp/2.0)+' nT', /norm

    endif
;------------------------------------------------------------------------

  endif

  if strlen(dmspfile) gt 0 then begin

    plot, [0,!d.x_size],[0,!d.y_size],/nodata,/noerase,		$
	xstyle=5, ystyle=5, pos = [0.0,0.0,1.0,1.0]

    loc = where(dmspdata(0,*) ne 0.0, count)
    if count gt 1 then begin

      fac = (float(!d.x_size)/10.0)/1000.0

      for i=loc(0)+1,n_elements(dmspdata(0,*))-1 do begin
        pointn=[dmspdata(0:1,i)]
        pointo=[dmspdata(0:1,loc(0))]
        length = fac*dmspdata(2,i)
        plot_dmsp_vector, pointn, pointo, length, /white
      endfor

    endif

    pointn=[float(!d.x_size)/100.0+1.0,float(!d.y_size)/100.0]
    pointo=[float(!d.x_size)/100.0,float(!d.y_size)/100.0]
    length = 1000.0*fac
    plot_dmsp_vector, pointn, pointo, length

    xyouts,0.02, 0.01, '1 Km/s',/norm

  endif

  plot, [-range, range], [-range, range],         $
        xstyle = 5, ystyle = 5,                 $
        /noerase, pos = cgmpos, /nodata

  if cgm then add_cgm, minlat, maxlat, rotation, pos

;  if r eq 0 then begin
    time_stamp = strmid(time,7,4)+':'+strmid(time,11,2)+' UT'
;  endif else begin
;    time_stamp = strmid(time,7,4)+':'+strmid(time,11,2)+' UT'+		$
;                 '  Residual'
;  endelse
  if rough eq 0 then 							$
    xyouts, range, range+range/100, time_stamp, alignment = 1.0
  if r ne 0.0 then $
;  xyouts, range, range-range/4, 'Residual', alignment = 1.0

;  xyouts, range, range+range/25, time_stamp, alignment = 1.0

  if (plot_num eq 1) and (rough eq 0) then color_bar, ncolors, maxval, r

  if (!d.name ne 'PS') and 						$
     (!d.name ne 'Z') and 						$
     (nn ne nfiles-1) then prompt_for_next

  if ((plot_num eq ppp-1) or 						$
     (nn eq nfiles-1)) and 						$
     (!d.name eq 'PS') then begin
    device, /close
    set_plot,'X'
  endif

  if !d.name eq 'Z' then begin
    bytemap = tvrd()
    tvlct, rr, gg, bb, /get
    if (plotloc eq 3) then begin
      if r eq 0 then giffile = time+'.gif' else giffile = time+'_sub.gif'
      write_gif, giffile, bytemap, rr, gg, bb
    endif else begin
      pic_image(plot_num,*,*) = bytemap
      time_st(plot_num) = time_stamp
      uttimesave(plot_num) = uttime
      if (plotloc eq 2) and 						$
	((plot_num eq ppp-1) or (nn eq nfiles-1)) and 			$
	(rough eq 1) then begin

        if r eq 0 then psfile = strmid(time,0,13)+'.ps'			$
        else psfile = strmid(time,0,13)+'_sub.ps'

        setdevice, psfile, 'l',4
        space = 0.04
        pos_space,ppp,space,sizes

        get_position,ppp,space,sizes,0,pos
        mvpos = pos(0)

	tp = min([nfiles-nn, ppp])

	for i=0,ppp-1 do begin
          get_position,ppp,space,sizes,i,pos
	  pos(0) = pos(0) - mvpos
	  pos(2) = pos(2) - mvpos
	  x_s = pos(2) - pos(0)
	  y_s = pos(3) - pos(1)
          !p.region = pos
          map_set,lat, long, rotation, limit = limit, 			$
	    /stereo, /noborder, /advance
	  loc = where(pic_image(i,*,*) gt 0,count)
	  if count gt 0 then 						$
	    tv, pic_image(i,*,*),pos(0),pos(1),xsize=x_s,ysize=y_s,/norm
          if (strmid(geolats,0,1) eq 'y') then begin
            MAP_CONTINENTS,MLINETHICK=1, color = 255
            np = maxlon - minlon + 1
            la = [fltarr(np)+minlat,fltarr(np)+maxlat,minlat]
            lo = [findgen(np)+minlon,maxlon-findgen(np),minlon]
            for ij=minlat,maxlat,10 do oplot, lo,la+ij-minlat, 		$
		linestyle=1, color=255
            la = [minlat,maxlat]
            lo = [0.0,0.0]
            for ij=0,360,45 do oplot, lo+ij,la,linestyle=1, color = 255
          endif
  	  if (strlen(magfile) gt 0) and (count gt 0) then begin

	    read_mag,magfile,uttimesave(i),magdata,names

; this is the geographic location of the magnetic pole for 1996 :
	    glat = 81.30
	    glon = 277.64	
	    if minlon lt 0.0 then glon = glon - 360
	    pole = convert_coord(glon, glat, /to_device)

	    for ij=0,n_elements(magdata(0,*))-1 do begin
	      if minlon lt 0.0 then magdata(1,ij)=magdata(1,ij)-360.0
   
	      if (magdata(0,ij) ge minlat) and (magdata(0,ij) le maxlat) and $
	         (magdata(1,ij) ge minlon) and (magdata(1,ij) le maxlon) then $
                dumb = convert_coord(magdata(1,ij), magdata(0,ij), /to_device)$
              else begin
                magdata(2,ij) = 0.0
                dumb = [0.0,0.0]
              endelse
              magdata(0,ij) = dumb(0)
              magdata(1,ij) = dumb(1)
            endfor

            plot, [0,!d.x_size],[0,!d.y_size],/nodata,/noerase,		$
	      xstyle=5, ystyle=5, pos = [0.0,0.0,1.0,1.0]

            loc = where((magdata(0,*) ne 0.0) and 			$
			(abs(magdata(2,*)) lt 2000.0), count)

            if count gt 1 then begin

              maxp=200.0	
              fac = abs((float(!d.x_size)/5.0)/maxp)

              for ij=loc(0),n_elements(magdata(0,*))-1 do begin

	        if magdata(2,ij) ne -1.0e32 then begin
	          symbol=names(ij)	
	          magpos = [magdata(0:1,ij)]
	          north = magdata(3,ij)
	          east = -1.0*magdata(2,ij)
                  plot_polar_vector, magpos, north, east, fac, pole, /white
                endif

              endfor

            endif else begin
              fac = abs((float(!d.x_size)/5.0)/2000.0)
              maxp = 200.0
            endelse

            pointo=[float(!d.x_size)/5.0,float(!d.y_size)/100.0+1.0]
            pointn=[float(!d.x_size)/5.0,float(!d.y_size)/100.0]

            length = maxp*fac/2.0
            plot_dmsp_vector, pointn, pointo, length

            xyouts,0.21, 0.02,  $ 
              '      ' +tostr(maxp/2.0)+' nT', /norm

          endif
;------------------------------------------------------------------------

	  px = pos(2)-0.01
          py = pos(3)+0.005
          xyouts, px, py, time_st(i), alignment = 1.0, 			$
		/norm, charsize = 0.9

          if (i eq 1) then 						$
	    plotct, ncolors, [0.97,0.3,1.0,0.7], maxval/10.0, 		$
		    'Photons/(cm2s)',/right

        endfor

      endif
    endelse

    device, /close
    set_plot,'X'
  endif

endfor

return

end

;-----------------------------------------------------------------------------
; main code
;-----------------------------------------------------------------------------

done = 0
time = 0

while (done eq 0) do begin

  if time eq 0 then begin
    polar, latitude, longitude, image, times
  endif else begin
    polar, latitude, longitude, image, times, /replot
  endelse

  print, '1. Use same data again for more plotting.'
  print, '2. Get new data set and plot.'
  print, '3. Exit program'
  print, 'Enter selection:'
  que=0
  read, que
  if que eq 1 then time = 1						$
  else if que eq 2 then delvar, latitude, longitude, image, times	$
  else done = 1

endwhile

end
