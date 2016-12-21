cur_t = -20
dt = 5
remote_dir = '/cloud/amie/html'

openw,2,'future_pot.html'
printf,2,'<title>'
printf,2,' '
printf,2,'The Future is HERE!'
printf,2,' '
printf,2,'</title>'
printf,2,' '
printf,2,'<body'
printf,2,' '
printf,2,'  bgcolor       = "#ffffff" '
printf,2,'  text          = "#000000" '
printf,2,'  link          = "#0000ff" '
printf,2,'  vlink         = "#ff9999"'
printf,2,' '
printf,2,'>'
printf,2,' '
printf,2,'<META HTTP-EQUIV="refresh" CONTENT="300">'
printf,2,' '
printf,2,' '
printf,2,'<center>'
printf,2,'<h1>'
printf,2,' '

openw,4,'ae.html'
printf,4,'<title>'
printf,4,' '
printf,4,'AE and Dst Estimates from AMIE'
printf,4,' '
printf,4,'</title>'
printf,4,' '
printf,4,'<body'
printf,4,' '
printf,4,'  bgcolor       = "#ffffff" '
printf,4,'  text          = "#000000" '
printf,4,'  link          = "#0000ff" '
printf,4,'  vlink         = "#ff9999"'
printf,4,' '
printf,4,'>'
printf,4,' '
printf,4,'<META HTTP-EQUIV="refresh" CONTENT="300">'
printf,4,' '
printf,4,' '
printf,4,'<center>'
printf,4,'<h1>'
printf,4,' '

openw,3,'rcp_script'
printf,3,'#!/bin/sh'
printf,3,'cp rtamie.html '+remote_dir
printf,3,'cp future_pot.html '+remote_dir
printf,3,'cp ae.html '+remote_dir

amie_file = 'output.binary'

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp

count = n_elements(time)
potential = reform(data(*,0,*,*))
fac       = reform(data(*,7,*,*))
by = reform(imf(*,1))
bz = reform(imf(*,2))

;--------------------------------------------------------------
; figure out grid:
;--------------------------------------------------------------

nlats = n_elements(lats)
nmlts = n_elements(mlts)

lat2d = fltarr(nmlts,nlats)
lon2d = fltarr(nmlts,nlats)

for i=0,nlats-1 do lon2d(*,i) = mlts*!pi/12.0 - !pi/2.0
for i=0,nmlts-1 do lat2d(i,*) = lats

xpos = (90.0-lat2d)*cos(lon2d)
ypos = (90.0-lat2d)*sin(lon2d)

for i=0,count-1 do begin

  if (i eq 0) then begin
    x_size = 550
    y_size = 450
  endif else begin
    x_size = 360
    y_size = 280
  endelse

  set_plot, 'Z', /copy, /interpolate
  device, set_resolution=[x_size,y_size], z= 0
  dy = float(!d.y_ch_size)/float(!d.y_size)

  plotdumb

  ppp   = 1
  space = 0.04
  pos_space,ppp,space,sizes
  get_position,ppp,space,sizes,0,pos
  pos([0,2]) = pos([0,2]) - pos(0) + space

  ctname = getenv('RTAMIE')+'/lib/idl/'+'blue_white_red.ct'
  readct,ncolors,ctname

  polyfill,[0.0,0.0,1.1,1.1,0.0],[0.0,1.1,1.1,0.0,0.0],color=ncolors

  mr = 40.0
  loc = where((90.0-lat2d(0,*)) le mr)

  maxi = 1.0e-6
  scale = maxi/14.0
  mini = -1.0*maxi

  levels = findgen(29)*scale - 14.0*scale
  clevels = fix(float(ncolors)*(findgen(29)+0.5)/30.0)
  facmini = min(levels)
  facmaxi = max(levels)

  contour, fac(i,*,loc), xpos(*,loc), ypos(*,loc), 		$
	nlevels=30, /noerase, pos = pos, xstyle = 5,            	$
	ystyle = 5, /cell_fill, levels=levels,				$
	c_colors = clevels,xrange = [-mr,mr],yrange = [-mr,mr]

  maxi = 75.0
  scale = float(round(maxi/14.0+0.5))
  mini = -1.0*scale*(float(round((maxi/scale)))-1.0)

  levels = findgen(29)*scale - 14.0*scale
  clevels = fix(float(ncolors)*(findgen(29)+0.5)/30.0)
  mini = min(levels)
  maxi = max(levels)

  contour, potential(i,*,loc), xpos(*,loc), ypos(*,loc), 		$
	nlevels=30, /noerase, pos = pos, xstyle = 5,            	$
	ystyle = 5, levels=levels, c_linestyle = 3.0*(levels lt 0.0),	$
	xrange = [-mr,mr],yrange = [-mr,mr], color = 0,/follow

  c_r_to_a, itime, time(i)
  h = float(itime(3))
  m = float(itime(4))

  rot = 90.0 - (h + m/60.0)*15.0

  ppp   = 1
  space = 0.0
  pos_space,ppp,space,sizes
  get_position,ppp,space,sizes,0,pos
  pos([0,2]) = pos([0,2]) - 0.068
  pos([1,3]) = pos([1,3]) + 0.01
  pos(2) = pos(2) - 0.03

  !p.region=pos
  map_set, 59.0,-90.0,rot,/ortho, /cont,/noerase, 			$
	limit=[50,-180,90,180], color = 5, /noborder

; plot station locations:

  if (i eq 0) then begin

    master='master_psi.dat'
    openr,1,master
    line = ''
    readf,1,line
    nstat = fix(strmid(line,0,3))

    readf,1,line

    lat = fltarr(nstat)
    lon = fltarr(nstat)
    list = intarr(nstat)
    slist = strarr(nstat)
    n = 0

    for j=0,nstat-1 do begin

      readf,1,line
      lat(j) = float(strmid(line,21,6))
      lon(j) = float(strmid(line,28,6))
      if (lat(j) gt 0.0) and (lat(j) lt 91.0) then begin
        list(n) = i
        slist(n) = strmid(line,4,2)
        n = n+1
      endif

    endfor

    close,1

    xyouts, lon, lat, slist, alignment=0.5, charsize = 0.7, color = 0

  endif

  ppp   = 1
  space = 0.04
  pos_space,ppp,space,sizes
  get_position,ppp,space,sizes,0,pos
  pos([0,2]) = pos([0,2]) - pos(0) + space

  plot, [-mr,mr],[-mr,mr],xstyle=5,ystyle = 5,/noerase,pos=pos,/nodata

  lats = 90.0-[fltarr(91)+(90.0-mr)-40.0, fltarr(91)+(90.0-mr), (90.0-mr)-30.0]
  lons = [findgen(91)*4.0,360.0-findgen(91)*4.0,0.0]*!pi/180.0

  x = lats*cos(lons)
  y = lats*sin(lons)
  polyfill, x,y, color = 101

  oplot, [-mr,mr],[mr,mr],color = 101, thick = 3
  oplot, [-mr,mr],[-mr,-mr],color = 101, thick = 3

  oplot, [mr,mr],[-mr,mr],color = 101, thick = 3
  oplot, [-mr,-mr],[-mr,mr],color = 101, thick = 3

  plotmlt, mr, /black

  xyouts, pos(0), pos(1), 'Min : '+tostr(round(min(potential(i,*,*))))+	$
	' kV', /norm, charsize = 0.8, color = 0
  xyouts, pos(0), pos(1)+dy, 'Max : '+tostr(round(max(potential(i,*,*))))+$
	' kV', /norm, charsize = 0.8, color = 0

  xyouts, pos(0), pos(3)-dy, 'By:'+tostr(fix(by(i))), /norm, charsize = 0.8, color = 0
  xyouts, pos(0), pos(3)-2.0*dy, 'Bz:'+tostr(fix(bz(i))), /norm, charsize = 0.8, color = 0

  c_a_to_s, itime, date
  time_str = strmid(date,10,5)
  date_str = strmid(date,0,9)
  xyouts, pos(2), pos(3)-dy, date_str, /norm, charsize = 0.8, 		$
	alignment = 1.0, color = 0
  xyouts, pos(2), pos(3)-2.0*dy, time_str, /norm, charsize = 0.8, 	$
	alignment = 1.0, color = 0

  xyouts, pos(2), pos(1), 'Contour:'+tostr(fix(scale))+' kV',	 	$
	/norm, charsize = 0.8, color = 0, alignment = 1.0

;  xyouts, mean(pos([0,2])), pos(1)+2.5*dy,'Southwest Research',	$
;	/norm,charsize=0.7,color=0,alignment=0.5

;  xyouts, mean(pos([0,2])), pos(1)+1.5*dy,'Institute',		$
;	/norm,charsize=0.7,color=0,alignment=0.5

  xyouts, mean(pos([0,2])), pos(1)+0.5*dy,'Aaron J. Ridley',		$
	/norm,charsize=0.7,color=0,alignment=0.5

  posct = [pos(2)+0.03,pos(1),pos(2)+0.06,pos(3)]
  plotct, ncolors, posct, [facmini,facmaxi], 'FAC (A/m)', /right, color = 0

  t = chopr('0'+tostr(itime(3)),2)+chopr('0'+tostr(itime(4)),2)
  if i eq 0 then t1 = t
  d = chopr('0'+tostr(itime(0)),2) +					$
      chopr('0'+tostr(itime(1)),2) +					$
      chopr('0'+tostr(itime(2)),2)

  giffile = 'amie.'+d+'.'+t1+'.'+t+'.gif'

  if i eq 0 then giffile_save = giffile

  b = tvrd()
  tvlct,rr,gg,bb,/get
  write_gif, giffile, b(*,0:n_elements(b(0,*))-2),rr,gg,bb

  if cur_t lt 0 then str = tostr(abs(cur_t))+' minutes ago<p>'		$
  else str = tostr(abs(cur_t))+' minutes from now<p>'
  cur_t = cur_t + dt
  if (i gt 0) then begin
    printf,2,'AMIE potential pattern for '+str
    printf,2,'<img src="'+giffile+'"><p>'
  endif
  printf,3,'cp '+giffile+' '+remote_dir
  printf,3,'/bin/rm '+giffile

  device, /close
  set_plot, 'X'

endfor

printf,2,'<CENTER><H1><FONT FACE="Arial, Helvetica">'
printf,2,'<a href="index.html"> Back To Main Page </a>'
printf,2,'</H1>'

dst = fltarr(1440)
ae = fltarr(1440)
t = findgen(1440)*60.0
dst_p = fltarr(1440)
stime = intarr(6)
etime = intarr(6)

openr,1,'ae_rt.final'

line = ''
readf,1,line
readf,1,line
readf,1,line
readf,1,line

done = 0
i = 0

while (not done) do begin

  readf,1,line
  if (eof(1)) then done = 1 else begin

    if (i eq 0) then for k=0,4 do stime(k) = fix(strmid(line,k*3,3))

    ae(i) = float(strmid(line,22,6))
    dst(i) = float(strmid(line,62,6))
    dst_p(i) = float(strmid(line,72,6))

    i = i + 1
    if (i eq 1440) then done = 1

  endelse

endwhile

close,1

c_a_to_r, stime, sut
eut = sut + 24.0*3600.0 - 60.0
c_r_to_a, etime, eut

mmdst = [-150.0,50.0]

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

set_plot, 'Z', /copy, /interpolate
device, set_resolution=[600,500], z= 0

ctname = getenv('RTAMIE')+'/lib/idl/'+'blue_white_red.ct'
readct,ncolors,ctname

plotdumb

polyfill,[0.0,0.0,1.1,1.1,0.0],[0.0,1.1,1.1,0.0,0.0],color=102

pos = [0.1,0.1,0.95,0.95]

plot, mm(t),[0,1],/nodata,xstyle=5,ystyle=5,pos=pos,/noerase

can_stime = 04.0*3600.0
can_etime = 13.0*3600.0
polyfill,[can_stime,can_stime,can_etime,can_etime,can_stime],		$
     [0,1,1,0,0], color = 90

eur_stime = 18.0*3600.0
eur_etime = 24.0*3600.0
polyfill,[eur_stime,eur_stime,eur_etime,eur_etime,eur_stime],		$
     [0,1,1,0,0], color = 90

loc = where(ae lt 5000,count)
if count eq 0 then begin
  ae = t*0.0
  loc = where(ae lt 5000,count)
endif

maxae = max(ae(loc))
maxi  = max([maxae,1000.0])

plot, t(loc), ae(loc), xstyle=1,			$
	ytitle = 'AE',				$
	xtickname = xtickname,			$
	xtitle = xtitle,			$
	xtickv = xtickv,			$
	xminor = xminor,			$
	xticks = xtickn, 			$
        color = 2,				$
	pos = pos, /noerase, yrange = [0,maxi],	$
	xrange = mm(t)

giffile = 'ae.'+d+'.'+t1+'.gif'

b = tvrd()
tvlct,rr,gg,bb,/get
write_gif, giffile, b(*,0:n_elements(b(0,*))-2),rr,gg,bb
printf,4,'AMIE Estimated AE for time '+t1+' UT<br>'
printf,4,'<img src="'+giffile+'"><br>'
printf,4,'Keep in mind that this only includes stations shown on the'
printf,4,'previous page.  If there are no stations in the midnight'
printf,4,'sectior, then the estimate is wrong!  Right now, we have'
printf,4,'good coverage in the midnight sector in the highlighted'
printf,4,'area of the plot.<p>'
printf,3,'cp '+giffile+' '+remote_dir
printf,3,'/bin/rm '+giffile

plotdumb
polyfill,[0.0,0.0,1.1,1.1,0.0],[0.0,1.1,1.1,0.0,0.0],color=50

plot, t, dst, xstyle=1, ystyle=1,		$
	yrange = mmdst, ytitle = 'DST',		$
	xtickname = xtickname,			$
	xtitle = xtitle,			$
	xtickv = xtickv,			$
	xminor = xminor,			$
	xticks = xtickn, color = 2, /noerase,	$
        thick = 2, pos = pos

oplot, t,dst_p, linestyle=2
oplot, [btr,etr],[0.0,0.0], linestyle = 1

giffile = 'dst.'+d+'.'+t1+'.gif'

b = tvrd()
tvlct,rr,gg,bb,/get
write_gif, giffile, b(*,0:n_elements(b(0,*))-2),rr,gg,bb
printf,4,'AMIE Estimated Dst for time '+t1+' UT<br>'
printf,4,'<img src="'+giffile+'"><br>'
printf,4,'This is 0 because we have NO low latitude stations at this time!<p>'
printf,3,'cp '+giffile+' '+remote_dir
printf,3,'/bin/rm '+giffile

printf,4,'<CENTER><H1><FONT FACE="Arial, Helvetica">'
printf,4,'<a href="index.html"> Back To Main Page </a>'
printf,4,'</H1>'

close,2,3,4
spawn,'chmod a+x rcp_script'

openw, 1, 'rtamie.html'

printf,1,'<title> '
printf,1,''
printf,1,'Real Time AMIE Page'
printf,1,''
printf,1,'</title>'
printf,1,''
printf,1,'<body'
printf,1,' '
printf,1,'  bgcolor	= "#ffffff" '
printf,1,'  text		= "#000000" '
printf,1,'  link		= "#0000ff" '
printf,1,'  vlink		= "#ff9999"'
printf,1,''
printf,1,'>'
printf,1,''
printf,1,'<META HTTP-EQUIV="refresh" CONTENT="300">'
printf,1,''
printf,1,''
printf,1,'<center>'
printf,1,'<h2>'
printf,1,''
printf,1,'rtAMIE potential pattern for 20 minutes ago<p>'
printf,1,''
printf,1,'<img src="'+giffile_save+'"> <p>'
printf,1,''
printf,1,'<a href="future_pot.html">Preditions for the next 60 minutes are HERE!</a><p>'
printf,1,'<a href="ae.html">rtAMIE derived AE and Dst are here.</a><p>'
printf,1,'</h2>'
printf,1,''
printf,1,'</center>'
printf,1,'The pattern above below the ionospheric potential pattern for 20 minutes'
printf,1,'ago.  It is based upon the IMF measured by the ACE'
printf,1,'satellite 80 minutes ago (60 minutes is the approximate time delay for the'
printf,1,'solar wind to propagate from the L1 point to the Earth), and a number of'
printf,1,'magnetometer measurements of the ionospheric current systems.  '
printf,1,'The magnetometer station locations are shown on the potential plot. <p>'
printf,1,''
printf,1,'Using the assimilative mapping of ionospheric electrodynamics'
printf,1,'(AMIE) (Richmond et al., JGR, 1988) technique to determine a large number'
printf,1,'of electric potential patterns from ground based magnetometer'
printf,1,'stations, and estimates of the hemispheric power index,  a'
printf,1,'multivariable linear regression analysis (e.g. Papitashvili et al.,'
printf,1,'1994) was preformed at each grid point to determine the relationship'
printf,1,'between the electric potential and the interplanetary magnetic field'
printf,1,'(IMF) Y and Z components.<p>'
printf,1,''
printf,1,'We have incorporated this statistically based model into the AMIE'
printf,1,'technique, allowing the statistical AMIE (SAMIE) to specify the'
printf,1,'background potential pattern, based upon measurements of the IMF.<p>'
printf,1,''
printf,1,'<CENTER><H1><FONT FACE="Arial, Helvetica">'
printf,1,'<a href=".."> Back To Main Page </a>'
printf,1,'</H1>'
printf,1,'</CENTER>'

close,1



end