pro plot_data_loc, namiein, date, time, dt, minlat

  scale_nt = 10.0
  scale_mv = 3.0

  print, date,' ',time

  mr = 90.0 - minlat

  sdate = strmid(date,4,2)+'-'+strmid(date,0,3)+'-'+strmid(date,10,2)+' '
  sdate = sdate + strmid(time,0,5)

  c_s_to_a, itime, sdate
  c_a_to_r, itime, rtime

  stime = rtime - dt
  etime = stime + 2.0*dt

  openr,11,namiein

  line = ''
  for i=1,4 do readf,11,line

  files =intarr(6)

  readf,11,files

  skip = files(0) + files(1) + files(2) + 1

  for i=1,skip do readf,11,line

  if (files(2) gt 0) then begin
    data_dir = line
    readf,11,line
  endif else readf,11,data_dir

  if (files(2) gt 0) then begin
    magfile = data_dir + strmid(line,8,16)
    magfile = strcompress(magfile,/remove_all)
  endif

  itime = intarr(6)
  latlon = fltarr(3)
  kloc = 0
  ang1n = 0.0
  en = 0.0
  ee = 0.0
  estring = ''

  t = findgen(10)*2.0*!pi/9
  r = 0.5
  usersym,r*cos(t),r*sin(t)

  psym = [8,4]

  for k=0,1 do begin

    for i=0,files(3+k)-1 do begin

      readf,11,line

      line = strcompress(line)
 
      for j=0,5 do line = strmid(line, strpos(line,' ')+1,strlen(line))
      filein = strmid(line,0,strpos(line,' '))

      print,"Searching file : ",filein

      filein = data_dir + filein

      openr,12,filein
      readf,12,line
      readf,12,line
      if (strpos(line,'year') lt 0) then readf,12,line

      done = 0

      while (not done) do begin
        readf,12,itime,latlon,kloc,ang1n,estring,			$
              format = '(I4,5I3,F6.2,F7.2,F6.2,I5,F7.2,A21)'
        c_a_to_r,itime, uttime
        if (uttime ge stime) and 					$
           (uttime le etime) and 					$
           (latlon(0) gt minlat) then begin
	  t = latlon(2)*2.0*!pi/24.0 - !pi/2
	  r = 90.0 - latlon(0)

          x0 = r*cos(t)
          y0 = r*sin(t)

          plots, [x0], [y0], psym = psym(k)

          ee = float(strmid(estring,0,7))
          ens = strcompress(strmid(estring,14,7))
          if strlen(ens) gt 1 then en = float(ens) else en = 0.0

          theta = (2.0*!pi - t) mod (2.0*!pi)
          d = (ee^2.0 + en^2.0)^0.5

; want pure north (h) to be 0 degrees, and rotate clock wise is positive

          phi = -1.0*asin(ee/d)
	  if (en lt 0.0) then phi = !pi - phi

          phi = phi - ang1n*!pi/180.0

          phi = !pi + (phi - theta)

          phi = phi + !pi

          xp = x0 + d*cos(phi)/scale_mv
          yp = y0 + d*sin(phi)/scale_mv

;if k eq 0 then $
;          plots, [x0,xp], [y0,yp]

        endif
        if (uttime gt etime) then done = 1
        if (eof(12)) then done = 1
      endwhile

      close,12

    endfor

  endfor

  done = 0
  cnt = 0

  while (not done) do begin
    readf,11,line
    if eof(11) then done = 1
    if (strpos(line, '&') ge 0) then cnt = cnt + 1
    if cnt eq 4 then done = 1
  endwhile

  t = findgen(10)*2.0*!pi/9
  r = 0.5
  usersym,r*cos(t),r*sin(t),/fill

  if (files(2) gt 0) then begin

    if (not eof(11)) then begin
      readf,11,line
      nstat = fix(strmid(line,0,3))
      readf,11,line
      lat = fltarr(nstat)
      lon = fltarr(nstat)
      weight = fltarr(nstat)
      sta = strarr(nstat)
      n = 0
      for i=0,nstat-1 do begin
        readf,11,line
        lat(n) = 90.0-float(strmid(line,36,6))
        lon(n) = float(strmid(line,42,7))*!pi/180.0
        weight(n) = float(strmid(line,68,4))
        sta(n) = strmid(line,4,3)
        n = n + 1 
      endfor

      lat = lat(0:n-1)
      lon = lon(0:n-1)
      sta = sta(0:n-1)

      c_r_to_a, itime, stime+dt
      h = float(itime(3))*60.0
      m = (float(itime(4))+h - 4.0*60.0 - 47.0)/4.0 * !pi/180.0

      lon = lon + m - !pi/2.0

      xpos = lat*cos(lon)
      ypos = lat*sin(lon)

      xpos2 = xpos
      ypos2 = ypos

      stime = stime + dt

      openr,12, magfile

      done = 0
      total = 0

      while (not done) do begin
        readf,12,line
        itime(0) = fix(strmid(line,0,2))  
        itime(1) = fix(strmid(line,2,2))
        itime(2) = fix(strmid(line,4,2))  
        itime(3) = fix(strmid(line,7,2))  
        itime(4) = fix(strmid(line,9,2))
        itime(5) = 0
        c_a_to_r,itime, uttime
        if (uttime ge stime) and 				$
           (uttime le etime) then begin
          etime = uttime
          stat_test = strmid(line,12,3)
          i = -1
          found = 0
          j = 0
          while (i eq -1) and (j lt n) do begin
            if (sta(j) eq stat_test) then begin
              i = j
              found = 1
            endif
            j = j + 1
          endwhile
          if (i eq -1) then done = 1				$
          else begin

	    total = total + 1
	    if total eq n then done = 1

	    h = float(strmid(line,19,6))
	    e = float(strmid(line,25,6))

            if (abs(h) lt 2000.0) and (abs(e) lt 2000.0) and	$
	       (lat(i) lt mr) and (weight(i) gt 0.0) then begin

              theta = (2.0*!pi - lon(i)) mod (2.0*!pi)
              d = (h^2.0 + e^2.0)^0.5

; want pure north (h) to be 0 degrees, and rotate clock wise is positive

              phi = -1.0*asin(e/d)
	      if (h lt 0.0) then phi = !pi - phi

              phi = !pi + (phi - theta)

; want to rotate counter clock wise 90 deg to the direction of convection

              phi = phi + !pi/2.0

              xpos2(i) = xpos(i) + d*cos(phi)/scale_nt
              ypos2(i) = ypos(i) + d*sin(phi)/scale_nt

              plots, [xpos(i)], [ypos(i)], psym = 8
              plots, [xpos(i),xpos2(i)], [ypos(i),ypos2(i)], thick=2

            endif

          endelse

        endif else for i=1,n-1 do readf,12,line

      endwhile

; plot scale

      xnt = 0.0
      ynt = -100.0

      x0 = mr
      y0 = mr

      plots, [x0], [y0], psym = 8
      plots, [x0,x0+xnt/scale_nt], [y0,y0+ynt/scale_nt]

      xyouts, mr-2.0, mr-2.0, '100 nT', alignment = 1.0, charsize = 0.6

    endif

    close,12

  endif

  close,11

  return

end


function sign, value

  loc = where(value ne 0.0,count)

  x = fltarr(n_elements(loc))+1.0

  if (count ne 0) then x(loc) = value(loc)/abs(value(loc)) 

  return,x

end

pro read_gang, filein, ltpos, lnpos, data, time, date, lats, lons, type

  type = 'Electric Potential'

  line = ''
  filein = filein+'.amiedat'

  year = 0
  print, 'Enter Year :'
  read, year
  if year lt 100 then year = year + 1900
  if year lt 1950 then year = year + 100
  year = tostr(year)

  month = 0
  print, 'Enter Month (1-12) :'
  read, month

  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'

  day = 0
  print, 'Enter Day :'
  read, day

  datei = strmid(mon,(month-1)*3,3)+' '+tostr(day)+', '+year

  close,1
  openr,1,filein

  readf,1, line
  readf,1, line
  lats = fix(strmid(line,6,3))

  readf,1, line
  readf,1, line
  lons = fix(strmid(line,7,3))
  for i=0,1 do readf,1,line
  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  for i=1,lons do for j=1,lats do begin
    ltpos(i-1,j-1) = float(j-1)*(90.0-50.0)/float(lats-1)
    lnpos(i-1,j-1) = 180.0 - float(i)*360.0/float(lons)
  endfor

  loc = where(lnpos lt 0.0, count)

  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  data = fltarr(1,lons,lats)
  fakedata = fltarr(1,6)
  dum = 1.0
  count = 0
  while (not eof(1)) do begin

    readf,1,line
    print, line
    hr = strmid(line,8,2)
    mi = float(strmid(line,10,3))*60.0
    readf,1,line
    if count gt 0 then begin
      data = [data, fltarr(1,lons,lats)]
      time = [time,hr+':'+chopr('0'+tostr(round(mi)),2)+' UT']
    endif else time = hr+':'+chopr('0'+tostr(mi),2)+' UT'

    for i=0,35 do begin
      for j=0,3 do begin
        readf,1, format = '(6E13.5)',fakedata
        data(count,i,j*6:(j+1)*6-1) = fakedata(0,*)
      endfor
      readf,1, format = '(1E13.5)',dum
      data(count,i,24) = dum
    endfor

    count = count + 1

  endwhile

  close,1

  date = strarr(count)

  for i=0,count-1 do date(i) = datei

  print, 'There were ',tostr(count),' times found'

  return

end

pro read_aaron, filein, ltpos, lnpos, data, time, date, lats, lons, type

  type = 'Electric Potential'

  line = ''
  fileti = filein+'.time'
  filein = filein+'.dat'

  close,1
  openr,1,filein

  readf,1, line
  readf,1, line
  lons = fix(strmid(line,13,3))

  for i=0,7 do readf,1,line
  readf,1, line
  lats = fix(strmid(line,12,3))
  for i=0,5 do readf,1,line
  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  for i=1,lons do for j=1,lats do begin
    ltpos(i-1,j-1) = float(j-1)*(90.0-50.0)/float(lats-1)
    lnpos(i-1,j-1) = 180.0-float(i)*360.0/float(lons)
  endfor

  loc = where(lnpos lt 0.0, count)

  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  data = fltarr(1,lons,lats)
  fakedata = fltarr(1,6)
  count = 0
  while (not eof(1)) do begin

    readf,1,line
    print, line
    if count gt 0 then data = [data, fltarr(1,lons,lats)]

    for i=0,24 do for j=1,36,6 do begin
      readf,1, format = '(6E13.5)',fakedata
      data(count,(j-1):(j+4),i) = fakedata(0,*)
    endfor

    count = count + 1

  endwhile

  close,1,2
  openr,2,fileti
  readf,2, line
  readf,2, line
  year = strmid(line,1,4)
  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'
  date = strarr(count)
  time = strarr(count)
  dt = fix(strmid(line,7,3))
  nt = fix(strmid(line,12,3))
  for i=0,count-1 do begin
    if (i eq 0 or nt eq 1) then begin
      readf,2,line
      month = fix(strmid(line,1,2))
      date(i) = strmid(mon,(month-1)*3,3)+' '+strmid(line,3,2)+', '+year
      h = fix(strmid(line,6,2))
      m = fix(strmid(line,8,2))
    endif else begin
      m = m + dt
      if (m ge 60) then begin
        m = m - 60
        h = h + 1
      endif
    endelse
    time(i) = chopr('0'+tostr(h),2)+':'+chopr('0'+tostr(m),2)+' UT'
  endfor
  close, 2

  return

end

pro read_ipot, filein, ltpos, lnpos, data, time, date, lats, lons, type

  type = 'Electric Potential'

  filelist = findfile(filein)
  count = n_elements(filelist)

  data = fltarr(count,25,34) 
  time = strarr(count)
  date = strarr(count)
  imf = strarr(count,3)

  for i=0,count-1 do begin

    data1 = read_izmem(filelist(i))
    data(i,0:23,1:33) = data1.pot
    data(i,24,1:33) = data(i,0,1:33)
    data(i,*,0) = mean(data(i,*,1))
    time(i) = data1.time + ' UT'
    imf(i,*) = data1.imf

    date(i) = data1.date
    itime = [fix(strmid(data1.date,2,2)),	$
	     fix(strmid(data1.date,5,2)),	$
	     fix(strmid(data1.date,8,2)),0,0,0]
    c_a_to_s, itime, sdate
    date(i) = strmid(sdate,3,3)+' '+		$
	      strmid(sdate,0,2)+', '+		$
	      strmid(data1.date,0,4)
  endfor

  lnpos = fltarr(25,34)
  ltpos = fltarr(25,34)
  for i=0,33 do lnpos(*,i) = findgen(25)*360.0/24.0 
  for i=0,24 do ltpos(i,*) = 90.0 - findgen(34)

  fileout = ''

  lats = 34
  lons = 25

  return

end

pro plotapot

  filein = ''
  fileout = ''
  line = ''

  read, 'Enter file name of the AMIE output : ', filein
  print, '1. AMIE output'
  print, '2. AMIE binary output'
  print, '3. IZMEM potential patterns'
  filetype = fix(ask('file type (1,2)','1'))

  que$ = ''
  read, 'Subtraction to be done (y/n) ? [n] ', que$
  if (strmid(que$,0,1) eq 'y') or (strmid(que$,0,1) eq 'Y') then begin
    sub = 1
    que$ = ''
    asub = 0
    print, '1. Subtract the previous pattern'
    print, '2. Subtract an averaged pattern'
    print, '3. Subtract a running average'
    read,'Enter option : ',asub
    asub = asub - 1
  endif else begin
    sub = 0
    asub = 0
  endelse
  fileout = ask('ps file name (return for screen)', '')

  ppp = 12
  ppp = fix(ask('number of plots per page',tostr(ppp)))

  pot = 5
  pot = float(ask('contour level spacing',tostr(pot)))

  if (filetype eq 1) then begin
    read_barbara, filein, ltpos, lnpos, data, time, date, lats, lons, 	$
	type, spd, by, bz
  endif

  if (filetype eq 2) then begin
    read_amie_binary, filein, data, lats, lons, time, type, imf,	$
                      ae, dst, hp, cpcp, date = date,			$
                      ltpos = ltpos, lnpos = lnpos, 			$
                      /plotapot, speed = spd, by = by, bz = bz
  endif

  if (filetype eq 3) then begin
    read_ipot, filein, ltpos, lnpos, data, time, date, lats, lons, type
  endif

  tmp_data = fltarr(lons,lats)

  maxran = 40.0

  que=''
  read, 'Would you like Images to be plotted under the convection? [n] ',que

  cplot = ask('whether you would like color or black and white (c or b)','b')

  if (mklower(que) eq 'y') then begin
    im = 1
    color = 255
    print, 'Enter the maximum value to the image (200 - 800) : '
    max_image = 0
    read, max_image
;------------------------------------------
    read_images, uttimes, image_save, time, max_image,ncolor
;------------------------------------------
  endif else begin
    im = 0
    if strlen(fileout) gt 0 then color = 0 else color = 255
  endelse

  que=''
  read, 'Would you like Data to be plotted under the convection? [n] ',que

  plot_data = 0
  if (mklower(que) eq 'y') then begin
    namiein = ask('namiein file (with full path)','./namiein.nov93')
    dt = float(ask('half window size (minutes)','05'))*60.0
    plot_data = 1
  endif

  if (sub eq 1) and (asub eq 1) then begin
    n1 = 0
    n2 = 0               
    print, 'Enter averaging start n, stop n :'
    read, n1,n2

    data(0,*,*) = data(n1,*,*)
    for n=n1+1,n2 do begin
      data(0,*,*) = data(0,*,*) + data(n,*,*)
    endfor
    data(0,*,*) = data(0,*,*)/(n2-n1+1)
  endif

  if (sub eq 1) and (asub eq 2) then begin
    runt = 0
    read, 'Enter time period of running average : ', runt
  endif

  print, 'Date of AMIE plots : ', date(0)
  print, 'Start time of file : ',time(0)
  print, 'End time of file : ',time(n_elements(time)-1)
  n2 = n_elements(time)
  print, 'Number of plots in File : ',n2
  n1 = fix(ask('starting plot number','0'))
  n2 = fix(ask('starting plot number',tostr(n2-1)))
  if n2 ge n_elements(time) then n2 = n_elements(time)-1

  if (n1 eq -1) then begin

    print, 'you have selected to enter individual numbers.'
    print, ''

    done = 0

    while (done eq 0) do begin

      que = ''
      read, 'Enter number (return to quit) : ',que
      if (strlen(que) eq 0) then done = 1			$
      else begin

        if (n_elements(list) eq 0) then list = fix(que)		$
	else list = [list,fix(que)]

      endelse

    endwhile

    n1 = 0
    n2 = n_elements(list)-1

  endif

  dmspfile = ''
  read, 'Enter the DMSP input file into AMIE (return for none): ', dmspfile

;------------------------------------------
  if strlen(dmspfile) gt 0 then begin
    if n_elements(list) gt 0 then n = list(0) else n = n1
    print, 'reading DMSP data from file : ',dmspfile
    readdmsp, dmspfile, date(0), time(n), maxran, dmspdata, /interpolate
    ppp = 2
  endif
;------------------------------------------

  if strlen(fileout) gt 0 then begin
     if ppp eq 9 then setdevice, fileout,'p',4,0.83		$
     else setdevice, fileout,'l',4,0.95
  endif else window,1,xsize=800,ysize=600

  plotdumb
  dy = float(!d.y_ch_size)/float(!d.y_size)

  space = 0.01
  space = space*9.0/float(ppp)
  pos_space, ppp, space, sizes

  aaa = !pi*findgen(361.0)/180.0
  xp = sin(aaa)
  yp = cos(aaa)
  maxran = 40.0

  post = fltarr(ppp,4)

  get_position,ppp,space,sizes,0,pos
  ys = 1.0-pos(3)
  for i=0,ppp-1 do begin
    get_position,ppp,space,sizes,i,pos
    if strlen(dmspfile) eq 0 then pos([1,3]) = pos([1,3]) + ys
    post(i,*) = pos
  endfor

  if strlen(fileout) eq 0 then potout = amiefile+'.pot'		$
  else potout = strmid(fileout,0,strpos(fileout,'.ps'))+'.pot'
  if (sub eq 1) then openw,1,potout

  lev = findgen(30)*pot
  if strpos(cplot,'c') gt -1 then begin
    levels = lev - 15.0*pot
    cfile = getenv('IDL_EXTRAS')+'blue_white_red.ct'
    readct, ncolors, cfile
    clevels = findgen(30)*ncolors/29.0
  endif

  n1sub = n1-sub
  if asub ne 1 then n1sub = n1
  if n1sub eq 0 and asub eq 0 and sub eq 1 then n1sub = 1

  xpos = (90.0-ltpos)*cos(lnpos*!pi/180.0 - !pi/2.0)
  ypos = (90.0-ltpos)*sin(lnpos*!pi/180.0 - !pi/2.0)
  dpos = (xpos^2.0 + ypos^2.0)^0.5

  for ij=n1sub,n2 do begin

    if (ij eq n1-1) then n = 0 else n = ij
    if asub eq 1 then subn = 0 else subn = ij-1

    if asub ne 1 then pn = (n-n1) mod ppp			$
    else if ij ge n1 then pn = (ij-(n1-sub)+2) mod ppp else pn = 1

    if strlen(dmspfile) gt 0 then pn = 0

    if (asub eq 2) then begin

      subn = 0
      ni = ij - runt/2
      if (ni le 0) then ni = 1
      nj = ni + runt
      if (nj gt n2) then nj = n2
      nc = nj - ni + 1

      data(0,*,*) = data(ni,*,*)
      for l=ni+1,nj do begin
        data(0,*,*) = data(0,*,*) + data(l,*,*)
      endfor
      data(0,*,*) = data(0,*,*)/float(nc)

    endif

    pos = post(pn,*)

;------------------------------------------
    if im eq 1 then begin
      tv, image_save(ij-n1sub,*,*), pos(0), pos(1), 		$
	xsize=pos(2)-pos(0), ysize=pos(3)-pos(1),/norm
    endif
;------------------------------------------

    if (n_elements(list) gt 0) and (ij ge 0) then n = list(ij)

    if ij eq n1-1 then begin

      plot, [0,1], /nodata, /noerase, pos = [0,0,1,1], 		$
	xstyle=5, ystyle=5
      if ppp eq 9 then begin
        xyouts, pos(0)-0.04, pos(3)-0.015, '(A)', charsize = 1.2, color = color
        xyouts, pos(0) - (pos(2)-pos(0)+space), pos(1)-space/2.0+0.005,'(B)', $
	  charsize=1.2, color = color
        oplot, [pos(0)-(pos(2)-pos(0)+space)+0.05,		$
		pos(2)+(pos(2)-pos(0)+space)], 			$
	       [pos(1)-space/2.0,pos(1)-space/2.0], color = color
      endif
    endif

    plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase

    xyouts, 0.5, 1.03, date(0), alignment = 0.5, 		$
	charsize=1.2, /norm

    xyouts, 0.96, 1.03, type, alignment = 1.0, charsize = 1.2, /norm
    xyouts, -0.01, 0.03, filein, charsize = 0.75, /norm, orient = 90

    if ij ne n1-1 then begin
      xyouts, -1.0*maxran, 0.90*maxran, time(n), 		$
        alignment = 0.0, charsize=0.9, color = color
    endif else begin
      xyouts,-0.95*maxran, 0.90*maxran, 'Steady', 		$
	alignment = 0.0, charsize=0.9, color = color
    endelse

    if (sub eq 1) and (ij ne n1-1) then begin
      xyouts, 0.95*maxran, 0.90*maxran, 'Residual',  	$
	alignment = 1.0, charsize=0.7, color = color
    endif

    tp = sizes.nby*sizes.nbx
    if (n2-ij lt sizes.nbx) or (pn+1 gt tp-sizes.nbx) then       $
      xyouts, 0.0,maxran*(-1.1), '00', alignment=0.5, 		$
	charsize=0.8
    if (pn mod sizes.nbx eq sizes.nbx-1) or (n2-ij eq 0) then   $
      xyouts, maxran*(1.01), -0.05*maxran, '06', 		$
	charsize=0.8
    if (pn lt sizes.nbx) then                 $
      xyouts, 0.0,maxran*(1.02), '12', alignment=0.5, 		$
	charsize=0.8
    if (pn mod sizes.nbx eq 0) or (ij eq n1) then                $
      xyouts, maxran*(-1.01), -0.05*maxran, '18', 		$
	alignment=1.0, charsize=0.8

    for i=0,lons-1 do for j=0,lats-1 do begin
      if (sub eq 1) and (ij ne n1-1) then 				$
        tmp_data(*,*) = (data(n,*,*)-data(subn,*,* ))			$
      else								$
        tmp_data(*,*) = data(n,*,*)
    endfor

    if strpos(cplot,'b') gt -1 then 					$
      levels = lev + float(fix(min(tmp_data)/pot))*pot - pot

    loc = where(dpos(0,*) le maxran,count)

    if count gt 0 then begin

      if ((!d.name eq 'PS') and (im eq 1)) or (!d.name eq 'X') then 	$
        contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
          pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
          xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
          color = 255, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
          levels = levels, nlevels = n_elements(levels)-1		$
      else begin

	if strpos(cplot,'c') gt -1 then begin
          contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
            pos = pos, /noerase, xstyle = 5, ystyle = 5,		$
            xrange = [-maxran,maxran], yrange = [-maxran,maxran],	$
            c_colors = clevels,	/cell_fill,                             $
            levels = levels, nlevels = n_elements(levels)-1
        endif

        contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
          pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
          xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
          color = 0, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
          levels = levels, nlevels = n_elements(levels)-1

      endelse

      if (strpos(mklower(type),'electric') gt -1) then begin
        clip = float(!p.clip)
        yoff = 0.5*maxran*float(!d.y_ch_size)/(clip(3)-clip(1))
        locm = where(tmp_data eq min(tmp_data))
        xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '-', alignment=0.5
        locm = where(tmp_data eq max(tmp_data))
        xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '+', alignment=0.5
      endif

      if max(tmp_data(*,loc)) lt 100 then begin
        xyouts,maxran*(1.0),maxran*(-0.96), 				$
	    'Ma: '+tostr(max(tmp_data(*,loc))), alignment=1.0, 		$
		  charsize=0.85, color = color
      endif else begin
        xyouts,maxran*(1.0),maxran*(-0.96), 				$
	  'Ma:'+tostr(max(tmp_data(*,loc))), alignment=1.0, 		$
		charsize=0.85, color = color
      endelse

      if (strpos(mklower(type),'cond') eq -1) then begin
        if min(tmp_data(*,loc)) gt -100 then begin
          xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	    'Mi: '+tostr(min(tmp_data(*,loc))), alignment=0.0, 		$
		  charsize=0.85, color = color
        endif else begin
          xyouts,maxran*(-1.005),maxran*(-0.96),	 		$
	    'Mi:'+tostr(min(tmp_data(*,loc))), alignment=0.0, 		$
		charsize=0.85, color = color
        endelse
      endif

      oplot, maxran*xp, maxran*yp, color = color
      oplot, 30.0*xp, 30.0*yp,linestyle=1, color = color
      oplot, 20.0*xp, 20.0*yp,linestyle=1, color = color
      oplot, 10.0*xp, 10.0*yp,linestyle=1, color = color
      oplot, [0.0,0.0],[-maxran,maxran], color = color
      oplot, [-maxran,maxran], [0.0,0.0], color = color

;      xyouts, pos(2), pos(3)-dy, by(ij), /norm, charsize = 0.7, 	$
;	      color = color, alignment = 1
;      xyouts, pos(2), pos(3)-1.7*dy, bz(ij), /norm, charsize = 0.7, 	$
; 	      color = color, alignment = 1

      if (plot_data) then plot_data_loc, namiein, 		$
        date(n), time(n), dt, 90.0-maxran

      if ((sub eq 1) and (ij gt n1-sub)) then 				$
        printf,1,min(tmp_data(*,loc)), max(tmp_data(*,loc))

    endif

;    if (pn+1 eq ppp) or (ij eq n2) then begin
;      xyouts, 0.5, -0.05, figure, /norm, alignment = 0.5
;    endif

;-------------------------------------------
    if (strlen(dmspfile) gt 0) then begin

      oplot, dmspdata(3,*),dmspdata(4,*), psym = 4

      pn = 1
      pos = post(pn,*)
      pos(0) = pos(0) + 0.05
      yr = pos(3) - pos(1)
      pos(1) = pos(3) - yr/2.0

      min1 = float(fix(min(tmp_data(*,loc))/10.0)-1)*10.0
      min2 = float(fix(min(dmspdata(2,*))/10.0)-1)*10.0
      max1 = float(fix(max(tmp_data(*,loc))/10.0)+1)*10.0
      max2 = float(fix(max(dmspdata(2,*))/10.0)+1)*10.0
      minmax = mm([min1,max1,min2,max2])

      nele = n_elements(dmspdata(1,*))
      lon  = fltarr(nele)
      rang = fltarr(nele)
      pot  = fltarr(nele)
      xamie = fltarr(nele)
      pamie = fltarr(nele)
      damie = fltarr(nele)
      lon(0:nele-1)  = 180-dmspdata(1,0:nele-1)*360.0/24.0 ; -180.0
      pot(0:nele-1) = dmspdata(2,0:nele-1)
      loc = where(lon lt 0.0,count)
      if count gt 0 then lon(loc) = lon(loc)+360.0
      rang(0:nele-1) = 90.0-dmspdata(0,0:nele-1)

      for i=0,nele-1 do begin
        difflon = abs(condata(0,*)-lon(i))
        difflat = abs(condata(1,*)-rang(i))
        loc = where(difflon eq min(difflon) and difflat eq min(difflat),count)
	xamie(i) = condata(1,loc(0))*sin(condata(0,loc(0))*!pi/180.0)
        pamie(i) = condata(2,loc(0))
        damie(i) = abs(condata(2,loc(0))-pot(i))
      endfor

      xdmsp = dmspdata(3,*)
      plot, xdmsp, dmspdata(2,*), xstyle = 1, 		$
	xrange = [-maxran, maxran], /noerase, pos = pos,	$
	yrange = minmax, ystyle = 1, ytitle = 'Potential (kV)',	$
	xtickname = strarr(10)+' '
      oplot, xamie, pamie, linestyle = 2

      mmr = minmax(1) - minmax(0)

      oplot, [-maxran,maxran],[0.0,0.0], linestyle = 1
      oplot, [-maxran+0.2*maxran,-maxran+0.4*maxran],		$
	[0.1*mmr+minmax(0),0.1*mmr+minmax(0)]
      xyouts, -maxran+0.5*maxran, 0.1*mmr+minmax(0), 'DMSP'

      oplot, [0.2*maxran,0.4*maxran],				$
	[0.1*mmr+minmax(0),0.1*mmr+minmax(0)], linestyle = 2
      if (filetype ne 4) then 					$
        xyouts, 0.5*maxran, 0.1*mmr+minmax(0), 'AMIE'		$
      else xyouts, 0.5*maxran, 0.1*mmr+minmax(0), 'IZMEM'

      sdtime = fix(dmspdata(5,0))*100
      edtime = fix(dmspdata(5,nele-1))*100

      minu = fix(60.0*float(fix(dmspdata(5,0)*100)-sdtime)/100.0)
      sdtime = sdtime + minu

      minu = fix(60.0*float(fix(dmspdata(5,nele-1)*100)-edtime)/100.0)
      edtime = edtime + minu

      xyouts, -maxran+0.1*maxran, 0.9*mmr+minmax(0), 		$
	'DMSP Start : '+tostr(sdtime)+' UT'
      xyouts, maxran-0.1*maxran, 0.9*mmr+minmax(0), 		$
	'DMSP End :'+tostr(edtime)+' UT', alignment = 1

      pos(1) = pos(3) - yr
      pos(3) = pos(3) - yr/2.0

      maxdiff = float(fix(max(damie)/10.0)+1)*10.0

      plot, xdmsp, damie, xstyle = 1, 				$
	xrange = [-maxran, maxran], /noerase, pos = pos,	$
	yrange = [0,maxdiff], ystyle = 1, 			$
	ytitle = 'Difference (kV)',				$
	xtitle = 'Degrees along Dawn - Dusk Meridian'
      oplot, [-maxran,maxran],[0.0,0.0], linestyle = 1

    endif
;-------------------------------------------

    if (pn+1 eq ppp) and (ij ne n2) then begin
      if !d.name eq 'X' then prompt_for_next
      plot, [0,1], /nodata, xstyle=5, ystyle=5

    endif

;------------------------------------------
      if ((im eq 1) and ((pn+1 eq ppp) or (ij eq n2))) then begin
	pos = [0.95,pos(1),1.0,pos(3)]
        color_bar, ncolor, pos, max_image, 'Reyleighs'
      endif
;------------------------------------------

  endfor

  if sub eq 1 then close,1

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end

