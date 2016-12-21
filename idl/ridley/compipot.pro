function sign, value

  loc = where(value ne 0.0,count)

  x = fltarr(n_elements(loc))+1.0

  if (count ne 0) then x(loc) = value(loc)/abs(value(loc)) 

  return,x

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


pro read_barbara, filein, ltpos, lnpos, data, time, date, lats, lons, type

  close,1
  openr,1,filein
 
  line = ''
  readf,1, line
  start = strpos(line, "AMIE")-20
  lons = fix(strmid(line,start+16,3))+1
  lats = fix(strmid(line,start+12,3))+1
  count = fix(strmid(line,start,4))
 
  format = '('+tostr(lons)+'E12.4,A5)'

  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  lts = fltarr(lats)
  fakedata = fltarr(1,lons)
  date = strarr(count)
  time = strarr(count)
  data = fltarr(count,lons,lats)
  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'
 
  done = 0
  while not done do begin
 
    readf,1, line
    if strpos(line,"ntime iyear") gt 0 then done = 1
 
  endwhile
 
  type = ''

  for i=0,count-1 do begin
 
    year = strmid(line,5,4)
    month = fix(strmid(line,10,2))
    date(i) = strmid(mon,(month-1)*3,3)+' '+strmid(line,13,2)+', '+year
    h = fix(strmid(line,16,2))
    m = fix(strmid(line,19,2))
    time(i) = chopr('0'+tostr(h),2)+':'+chopr('0'+tostr(m),2)+' UT'
 
    readf,1, line
    readf,1, line

    for j=0,lats-1 do begin
 
      readf,1, line
      lts(j) = float(strmid(line,5,5))
      readf,1, format=format,fakedata, type
      data(i,0:lons-1,j) = fakedata(0,0:lons-1)
 
    endfor
 
    if i lt count-1 then readf,1, line
 
  endfor

  close,1

  for i=0,lons-1 do for j=0,lats-1 do begin
    ltpos(i,j) = 90.0 - lts(j)
    lnpos(i,j) = float(i)*360.0/float(lons-1)
  endfor
 
  loc = where(lnpos lt 0.0, count)
 
  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  type = strmid(type,1,strlen(type)-1)

  if (strpos(type,'pot') gt -1) or (strpos(type,'efpo') gt -1) then begin
    data = data/1000.0
    type = 'Electric Potential'
  endif
  if strpos(type,'sigp') gt -1 then type = 'Pedersen Conductance'
  if strpos(type,'sigh') gt -1 then type = 'Hall Conductance'
  if strpos(type,'dfac') gt -1 then type = 'Field Aligned Current'
  if strpos(type,'sjht') gt -1 then begin
    type = 'Simple Joule Heating'
    data = data*1000.0                           ; mW/m^2 * 1000
  endif

  return

end

pro compipot

  line = ''

  filein1 = ask('file name of the 1st AMIE output','x950324.0000.nc.pot')
  title1 = ask('title for this file','AMIE')
  read_barbara, filein1, ltpos, lnpos, data1, time1, date1, lats1, lons1,type

  data1 = data1

  filein2 = ask('file name of the IZMEM output',			$
  		'd.izmem/iz*.pot')
  title2 = ask('title for this file','IZMEM')
  read_ipot, filein2, ltposi, lnposi, data2, time2, date2, lats2, lons2, type

  datai = fltarr(n_elements(time1),lons1,lats1)
  loc = indgen(lats1-1)*2
  datai(*,*,indgen(lats1-1)) = data2(*,*,loc)
  datai(*,*,lats1-1) = 0.0

  lats2 = lats1
  lons2 = lons1
  data2 = datai

  if lats1 ne lats2 then print, "Number of latitudes are not equal!"
  if lons1 ne lons2 then print, "Number of longitudes are not equal!"

  lats = lats1
  lons = lons1
  date = date1
  time = time1

  n = n_elements(time)-1

  fileout = ask('ps file name (return for screen) ', '')

  ppp = 12

  pot = 2
  pot = float(ask('contour level spacing',tostr(pot)))

  tmp_data = fltarr(lons,lats)
  aerrdata = fltarr(lons,lats)
  merrdata = fltarr(lons,lats)
  adata1   = fltarr(lons,lats)
  adata2   = fltarr(lons,lats)

  maxran = 40.0
  color = 255
  if strlen(fileout) gt 0 then begin
     if ppp eq 9 then setdevice, fileout,'p',4,0.83		$
     else setdevice, fileout,'l',4,0.95
     color = 0
     fileout = strmid(fileout,0,strpos(fileout,'.ps'))+'.diff.ps'
  endif else window,1,xsize=800,ysize=600

  plotdumb

  xpos = (90.0-ltpos)*cos(lnpos*!pi/180.0 - !pi/2.0)
  ypos = (90.0-ltpos)*sin(lnpos*!pi/180.0 - !pi/2.0)
  dpos = (xpos^2.0 + ypos^2.0)^0.5

  print, 'Date of AMIE plots : ', date(0)
  print, 'Start time of file : ',time(0)
  print, 'End time of file : ',time(n_elements(time)-1)
  n2 = n_elements(time)
  print, 'Number of plots in File : ',n2
  n1 = fix(ask('starting plot number','0'))
  n2 = fix(ask('starting plot number',tostr(n2-1)))
  if n2 ge n_elements(time) then n2 = n_elements(time)-1

  npts = n2 - n1 + 1
  cpcp = fltarr(4,npts)  

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
    post(i,*) = pos
  endfor

  pot2 = 3.0*pot
  lev = findgen(30)*pot
  lev2 = findgen(30)*pot2

  for ij=n1,n2 do begin

    n = ij
    pn = (n-n1) mod ppp

    pos = post(pn,*)

    plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase

    oplot, maxran*xp, maxran*yp, color = color
    oplot, 30.0*xp, 30.0*yp,linestyle=1, color = color
    oplot, 20.0*xp, 20.0*yp,linestyle=1, color = color
    oplot, 10.0*xp, 10.0*yp,linestyle=1, color = color
    oplot, [0.0,0.0],[-maxran,maxran], color = color
    oplot, [-maxran,maxran], [0.0,0.0], color = color

    xyouts, 0.5, 1.03, date(0), alignment = 0.5, 			$
	charsize=1.2, /norm

    xyouts, 0.96, 1.03, type, alignment = 1.0, charsize = 1.2, /norm
    xyouts, -0.03, 0.03, filein1+' - '+filein2, 			$
      charsize = 0.75, /norm, orient = 90

    xyouts, -1.0*maxran, 0.90*maxran, time(n), 				$
        alignment = 0.0, charsize=0.9, color = color

    tp = sizes.nby*sizes.nbx
    if (n2-ij lt sizes.nbx) or (pn+1 gt tp-sizes.nbx) then       	$
      xyouts, 0.0,maxran*(-1.1), '00', alignment=0.5, 			$
	charsize=0.8
    if (pn mod sizes.nbx eq sizes.nbx-1) or (n2-ij eq 0) then   	$
      xyouts, maxran*(1.01), -0.05*maxran, '06', 			$
	charsize=0.8
    if (pn lt sizes.nbx) then                 				$
      xyouts, 0.0,maxran*(1.02), '12', alignment=0.5, 			$
	charsize=0.8
    if (pn mod sizes.nbx eq 0) or (ij eq n1) then                	$
      xyouts, maxran*(-1.01), -0.05*maxran, '18', 			$
	alignment=1.0, charsize=0.8

    tmp_data(*,*) = data1(n,*,*)-data2(n,*,*)
    aerrdata      = aerrdata + abs(tmp_data)
    merrdata      = merrdata + tmp_data
    adata1        = adata1   + data1(n,*,*)
    adata2        = adata2   + data2(n,*,*)

    cpcp(0,n-n1) = max(data1(n,*,*))-min(data1(n,*,*))
    cpcp(1,n-n1) = max(data2(n,*,*))-min(data2(n,*,*))
    cpcp(2,n-n1) = max(tmp_data)-min(tmp_data)
    cpcp(3,n-n1) = mean(tmp_data)

    levels = lev + float(fix(min(tmp_data)/pot))*pot - pot

    loc = where(dpos(0,*) le maxran,count)

    if count gt 0 then begin

      if (!d.name eq 'X') then 						$
        contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
          pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
          xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
          color = 255, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
          levels = levels, nlevels = n_elements(levels)-1		$
      else								$
        contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
          pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
          xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
          color = 0, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
          levels = levels, nlevels = n_elements(levels)-1

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

    endif

    if (pn+1 eq ppp) and (ij ne n2) then begin
      if !d.name eq 'X' then prompt_for_next
      plot, [0,1], /nodata, xstyle=5, ystyle=5
    endif

  endfor

  closedevice

  if strlen(fileout) gt 0 then setdevice, fileout,'l',4,0.95

  plotdumb

  sdate = strmid(date(n1),4,2)+'-'+mklower(strmid(date(n1),0,3))+'-'+	$
	  strmid(date(n1),10,2) + ' '+ strmid(time(n1),0,5)
  edate = strmid(date(n2),4,2)+'-'+mklower(strmid(date(n2),0,3))+'-'+	$
	  strmid(date(n2),10,2) + ' '+ strmid(time(n2),0,5)
  c_s_to_a, itime, sdate
  c_a_to_r, itime, stime
  c_s_to_a, itime, edate
  c_a_to_r, itime, etime

  time_axis, stime, etime, s_time_range, e_time_range,        		$
        xtickname, xtitle, xtickvalue, xminor, xtickn

  times = findgen(npts)*(e_time_range-s_time_range)/float(npts-1)

  poss = fltarr(4)

  get_position, ppp, space, sizes, 0, pos
  plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase, /nodata
  tmp_data = adata1/float(npts)
  levels = lev2 + float(fix(min(tmp_data)/pot2))*pot2 - pot2
  loc = where(dpos(0,*) le maxran,count)
  if count gt 0 then begin
    contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
      pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
      xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
      color = color, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
      levels = levels, nlevels = n_elements(levels)-1
    if (strpos(mklower(type),'electric') gt -1) then begin
      clip = float(!p.clip)
      yoff = 0.5*maxran*float(!d.y_ch_size)/(clip(3)-clip(1))
      locm = where(tmp_data eq min(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '-', alignment=0.5
      locm = where(tmp_data eq max(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '+', alignment=0.5
    endif
    xyouts,maxran*(1.0),maxran*(-0.96), 				$
      'Ma: '+tostr(max(tmp_data(*,loc))), alignment=1.0, 		$
      charsize=0.85, color = color
    if (strpos(mklower(type),'cond') eq -1) then begin
      xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	'Mi: '+tostr(min(tmp_data(*,loc))), alignment=0.0, 		$
        charsize=0.85, color = color
    endif
  endif
  plotmlt,maxran,/no00,/no06
  xyouts,maxran*(-1.0),maxran*0.96, title1, color = color
  poss(0) = pos(0)
  poss(3) = pos(1) - space
  poss(1) = poss(3) - (pos(3)-pos(1))

  get_position, ppp, space, sizes, 1, pos
  plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase, /nodata
  tmp_data = adata2/float(npts)
  levels = lev2 + float(fix(min(tmp_data)/pot2))*pot2 - pot2
  loc = where(dpos(0,*) le maxran,count)
  if count gt 0 then begin
    contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
      pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
      xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
      color = color, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
      levels = levels, nlevels = n_elements(levels)-1
    if (strpos(mklower(type),'electric') gt -1) then begin
      clip = float(!p.clip)
      yoff = 0.5*maxran*float(!d.y_ch_size)/(clip(3)-clip(1))
      locm = where(tmp_data eq min(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '-', alignment=0.5
      locm = where(tmp_data eq max(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '+', alignment=0.5
    endif
    xyouts,maxran*(1.0),maxran*(-0.96), 				$
      'Ma: '+tostr(max(tmp_data(*,loc))), alignment=1.0, 		$
      charsize=0.85, color = color
    if (strpos(mklower(type),'cond') eq -1) then begin
      xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	'Mi: '+tostr(min(tmp_data(*,loc))), alignment=0.0, 		$
        charsize=0.85, color = color
    endif
  endif
  plotmlt,maxran,/no06,/no00,/no18
  xyouts,maxran*(-1.0),maxran*0.96, title2, color = color

  get_position, ppp, space, sizes, 2, pos
  plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase, /nodata
  tmp_data = merrdata/float(npts)
  levels = lev + float(fix(min(tmp_data)/pot))*pot - pot
  if count gt 0 then begin
    contour, tmp_data(*,loc), xpos(*,loc), ypos(*,loc), 		$
      pos = pos, /noerase, xstyle = 5, ystyle = 5,			$
      xrange = [-maxran,maxran], yrange = [-maxran,maxran],		$
      color = color, c_linestyle = 3.0*(levels lt 0.0), /follow,	$
      levels = levels, nlevels = n_elements(levels)-1
    if (strpos(mklower(type),'electric') gt -1) then begin
      clip = float(!p.clip)
      yoff = 0.5*maxran*float(!d.y_ch_size)/(clip(3)-clip(1))
      locm = where(tmp_data eq min(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '-', alignment=0.5
      locm = where(tmp_data eq max(tmp_data))
      xyouts, xpos(locm(0)), ypos(locm(0))-yoff, '+', alignment=0.5
    endif
    xyouts,maxran*(1.0),maxran*(-0.96), 				$
      'Ma: '+tostr(max(tmp_data(*,loc))), alignment=1.0, 		$
      charsize=0.85, color = color
    if (strpos(mklower(type),'cond') eq -1) then begin
      xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	'Mi: '+tostr(min(tmp_data(*,loc))), alignment=0.0, 		$
        charsize=0.85, color = color
    endif
  endif
  poss(2) = pos(2)
  plotmlt,maxran,/no00,/no18
  xyouts,maxran*(1.0),maxran*0.96, 'Difference', color = color, alignment = 1.0

  if (strpos(type,'Electric') ne -1) then ytitle = 'Potential (kV)'
  if (strpos(type,'Hall') ne -1) then ytitle = 'Hall Cond. (mhos)'
  if (strpos(type,'Ped') ne -1) then ytitle = 'Ped. Cond. (mhos)'
  if (strpos(type,'Field') ne -1) then ytitle = 'FAC (mA/m!X2!N)'
  if (strpos(type,'Joule') ne -1) then ytitle = 'Joule Ht. (mW/m!X2!N)*1000'

  plot, times, cpcp(0,*), pos = poss, /noerase, 			$
	xtickname = xtickname, xstyle = 1, 				$
	yrange = mm(cpcp(0:2,*)), ytitle = ytitle,			$
	xtickv=xtickvalue, xticks = xtickn, xminor = xminor, 		$
	xtitle = xtitle
  oplot, times, cpcp(1,*), linestyle = 2
  oplot, times, cpcp(2,*), linestyle = 1
  posc = [poss(2)+0.01, (poss(1)+poss(3))/2.0, poss(2)+0.06, poss(3)]
  plot, [0,1], /nodata, /noerase, xstyle=5, ystyle=5, pos = posc
  oplot, [0,1], [0.5,0.5]
  xyouts, 1.03, 0.5, title1
  oplot, [0,1], [0.25,0.25], linestyle = 2
  xyouts, 1.03, 0.25, title2
  oplot, [0,1], [0.0,0.0], linestyle = 1
  xyouts, 1.03, 0.0, 'Difference'

  closedevice

  return

end

