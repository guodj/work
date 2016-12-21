file = ask('netcdf file to plot','test_output_1.nc')

filelist = findfile(file)

nfiles = n_elements(filelist)

nt = 0

nstart = intarr(nfiles)
nend = intarr(nfiles)

for i=0,nfiles-1 do begin

  file = filelist(i)

  id = ncdf_open(file)

  ncdf_varget, id, 36, potential

  nstart(i) = nt
  nt = nt + n_elements(potential(0,0,0,*))
  nend(i) = nt-1

  if i eq 0 then begin
    ncdf_varget, id, 0, lon_save
    ncdf_varget, id, 1, lat_save
  endif

  ncdf_close,id

endfor

nlat = n_elements(potential(0,*,0,0))
nlon = n_elements(potential(*,0,0,0))
nalt = n_elements(potential(0,0,*,0))

altitude  = fltarr(nlon,nlat,nalt,nt)
otwo      = fltarr(nlon,nlat,nalt,nt)
oxygen    = fltarr(nlon,nlat,nalt,nt)
nitox     = fltarr(nlon,nlat,nalt,nt)
potential = fltarr(nlon,nlat,nalt,nt)
electron  = fltarr(nlon,nlat,nalt,nt)
un        = fltarr(nlon,nlat,nalt,nt)
vn        = fltarr(nlon,nlat,nalt,nt)
n4s       = fltarr(nlon,nlat,nalt,nt)
n2d       = fltarr(nlon,nlat,nalt,nt)
te        = fltarr(nlon,nlat,nalt,nt)
tn        = fltarr(nlon,nlat,nalt,nt)
ut        = fltarr(nt)

print, 'I found ',tostr(fix(nt)),' times in the file(s)'


for i=0,nfiles-1 do begin

  file = filelist(i)
  print,'Reading file : ',file

  id = ncdf_open(file)
  ncdf_varget, id, 8, ut_
  ncdf_varget, id, 35, altitude_
;  ncdf_varget, id, 24, otwo_
;  ncdf_varget, id, 25, oxygen_
;  ncdf_varget, id, 27, nitox_
  ncdf_varget, id, 36, potential_
  ncdf_varget, id, 32, electron_
  ncdf_varget, id, 22, un_
  ncdf_varget, id, 23, vn_
;  ncdf_varget, id, 26, n4s_
;  ncdf_varget, id, 29, n2d_
;  ncdf_varget, id, 31, te_
  ncdf_varget, id, 21, tn_

  ncdf_close,id

  altitude(*,*,*,nstart(i):nend(i))  = altitude_
;  otwo(*,*,*,nstart(i):nend(i))      = otwo_
;  oxygen(*,*,*,nstart(i):nend(i))    = oxygen_
;  nitox(*,*,*,nstart(i):nend(i))     = nitox_
  potential(*,*,*,nstart(i):nend(i)) = potential_
  electron(*,*,*,nstart(i):nend(i))  = electron_
  un(*,*,*,nstart(i):nend(i))        = un_
  vn(*,*,*,nstart(i):nend(i))        = vn_
;  n4s(*,*,*,nstart(i):nend(i))       = n4s_
;  n2d(*,*,*,nstart(i):nend(i))       = n2d_
;  te(*,*,*,nstart(i):nend(i))        = te_
  tn(*,*,*,nstart(i):nend(i))        = tn_
  ut(nstart(i):nend(i)) = ut_

endfor

print, ut

alts = fltarr(nalt,nt)

for j=0,nt-1 do $
  for i=0,nalt-1 do $
    alts(i,j) = total(altitude(*,*,i,j))/n_elements(altitude(*,*,i,j))/100000.0

lats = fltarr(nlon+1,nlat)
lons = fltarr(nlon+1,nlat)
la = (findgen(nlat)*5.0 + 2.5)
lo = findgen(nlon+1)*5.0 * !pi/180.0 + !pi/2.0

for i=0,nlat-1 do lons(*,i) = lo
for i=0,nlon do lats(i,*) = la

x = lats*cos(lons)
y = lats*sin(lons)

mr = 40

loc = where(la lt mr)

;setdevice

closedevice
setdevice, 'wind_s.ps', 'l', 4, 0.9

ppp = 6
space = 0.02
pos_space, ppp,space, sizes

plotdumb

  a = 9

  get_position, ppp, space, sizes, 0, pos

  pos_save = fltarr(4)
  pos_save(0) = pos(0)

  vel_north = -reform(vn(*,loc,a,0))
  vel_east  = reform(un(*,loc,a,0))

  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /long

  get_position, ppp, space, sizes, 1, pos

  vel_north = -reform(vn(*,loc,a,nt/2))
  vel_east  = reform(un(*,loc,a,nt/2))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	title = 'Neutral Wind at '+tostr(fix(alts(a,nt/2)))+' Km', /noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /no18, /long

  get_position, ppp, space, sizes, 2, pos
  pos_save(2) = pos(2)

  vel_north = -reform(vn(*,loc,a,nt-1))
  vel_east  = reform(un(*,loc,a,nt-1))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no18, /no00, /long

  get_position, ppp, space, sizes, 3, pos
  pos_save(1) = pos(1)
  pos_save(3) = pos(3)

  pos = pos_save

  vtotal = sqrt(vn^2 + un^2)/100.0
  plot, ut*60.0, vtotal(56,1,a,*), ytitle = 'm/s', pos = pos, /noerase, $
	yrange=mm(vtotal(56,1,a,*))

closedevice


la = 180.0 - (findgen(nlat)*5.0 + 2.5)
lo = findgen(nlon+1)*5.0 * !pi/180.0 + !pi/2.0

for i=0,nlat-1 do lons(*,i) = lo
for i=0,nlon do lats(i,*) = la

x = lats*cos(lons)
y = lats*sin(lons)

mr = 40

loc = where(la lt mr)

closedevice
setdevice, 'wind_n_dn.ps', 'p', 4, 0.9

ppp = 6
space = 0.02
pos_space, ppp,space, sizes, ny = 2

plotdumb

  a = 9

  get_position, ppp, space, sizes, 0, pos

  pos_save = fltarr(4)
  pos_save(0) = pos(0)

  it = [12,25,31]

  vel_north = reform(vn(*,loc,a,it(0)))
  vel_east  = reform(un(*,loc,a,it(0)))

  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /long

  get_position, ppp, space, sizes, 1, pos

  vel_north = reform(vn(*,loc,a,it(1)))
  vel_east  = reform(un(*,loc,a,it(1)))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  xyouts, (pos(0)+pos(2))/2.0, pos(3) + 2*space, $
	'Neutral Wind at '+tostr(fix(alts(a,nt/2)))+' Km', charsize = 1.25,$
	/norm,alignment = 0.5

  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /no18, /long

  get_position, ppp, space, sizes, 2, pos
  pos_save(2) = pos(2)

  vel_north = reform(vn(*,loc,a,it(2)))
  vel_east  = reform(un(*,loc,a,it(2)))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no18, /no00, /long

  get_position, ppp, space, sizes, 3, pos
  pos_save(1) = pos(1)
  pos_save(3) = pos(3)

  pos = pos_save

  vtotal = sqrt(vn^2 + un^2)/100.0
  plot, ut,vtotal(18,29,a,*), ytitle = 'm/s', pos = pos, /noerase, $
	yrange=mm(vtotal(18,29,a,*)), xtitle = 'UT Hours'

  oplot, [ut(it(0)),ut(it(0))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(1)),ut(it(1))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(2)),ut(it(2))], [0.0,10000.0], thick = 0.1

  xyouts, 0.1,167.0,'Winds at '+tostr(360+lon_save(18))+'!Eo!N lon by '+$
          tostr(lat_save(29))+'!Eo!N lat', charsize = 0.75

closedevice

closedevice
setdevice, 'wind_n_up.ps', 'p', 4, 0.9

ppp = 6
space = 0.02
pos_space, ppp,space, sizes, ny = 2

plotdumb

  a = 15

  get_position, ppp, space, sizes, 0, pos

  pos_save = fltarr(4)
  pos_save(0) = pos(0)

  it = [12,25,31]

  vel_north = reform(vn(*,loc,a,it(0)))
  vel_east  = reform(un(*,loc,a,it(0)))

  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /long

  get_position, ppp, space, sizes, 1, pos

  vel_north = reform(vn(*,loc,a,it(1)))
  vel_east  = reform(un(*,loc,a,it(1)))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  xyouts, (pos(0)+pos(2))/2.0, pos(3) + 2*space, $
	'Neutral Wind at '+tostr(fix(alts(a,nt/2)))+' Km', charsize = 1.25,$
	/norm,alignment = 0.5

  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no06, /no00, /no18, /long

  get_position, ppp, space, sizes, 2, pos
  pos_save(2) = pos(2)

  vel_north = reform(vn(*,loc,a,it(2)))
  vel_east  = reform(un(*,loc,a,it(2)))

;  scale = 10.0/max([vel_north,vel_east])

  plot, [-mr,mr],[-mr,mr],xstyle = 5, ystyle = 5, pos = pos, /nodata, 	$
	/noerase
  plot_vectors_polar, vel_east, vel_north, $
                      lats(0:nlon-1,loc), lons(0:nlon-1,loc),$
	              scale, /arrow, length = 10000.0, $
                      pos = [mr,mr], title = '100 m/s'
;, nskip = 2
  plotmlt,mr, /no18, /no00, /long

  get_position, ppp, space, sizes, 3, pos
  pos_save(1) = pos(1)
  pos_save(3) = pos(3)

  pos = pos_save

  vtotal = sqrt(vn^2 + un^2)/100.0
  plot, ut,vtotal(18,32,a,*), ytitle = 'm/s', pos = pos, /noerase, $
	yrange=mm(vtotal(18,32,a,*)), xtitle = 'UT Hours'

  oplot, [ut(it(0)),ut(it(0))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(1)),ut(it(1))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(2)),ut(it(2))], [0.0,10000.0], thick = 0.1

  xyouts, 0.1,375.0,'Winds at '+tostr(360+lon_save(18))+'!Eo!N lon by '+$
          tostr(lat_save(32))+'!Eo!N lat', charsize = 0.75

closedevice

setdevice, 'potential.ps', 'p', 4, 0.90

ppp = 6
space = 0.02
pos_space, ppp,space, sizes, ny = 2

readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
c_colors = (ncolors-1)*findgen(30)/29.0 + 1

plotdumb

ialt = 5

get_position, ppp, space, sizes, 0, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pos_save = fltarr(4)
pos_save(0) = pos(0)
pos_save(2) = pos(2)

maxi = max(potential(*,*,ialt,*))/1000.0
mini = min(potential(*,*,ialt,*))/1000.0
maxi = max([maxi,-mini])
mini = -maxi

levels = (maxi-mini)*findgen(9)/8.0 + mini
c_levels = (maxi-mini)*findgen(30)/29.0 + mini

it = [12,25,31]

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = potential(*,*,ialt,it(0))/1000.0
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /long

get_position, ppp, space, sizes, 1, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = potential(*,*,ialt,it(1))/1000.0
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long


xyouts, (pos(0)+pos(2))/2.0, pos(3) + 2*space, $
	'Potential at '+tostr(fix(alts(a,nt/2)))+' Km', charsize = 1.25,$
	/norm,alignment = 0.5

get_position, ppp, space, sizes, 2, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = potential(*,*,ialt,it(2))/1000.0
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'kV', /right, color = 0

get_position, ppp, space, sizes, 3, pos
pos_save(0) = pos(0) - space/2.0

get_position, ppp, space, sizes, 5, pos
pos_save(2) = pos(2) - space/2.0

pos_save(1) = pos(1)
pos_save(3) = pos(3)

pos = pos_save

cpcp = fltarr(nt)

for i=0,nt-1 do cpcp(i) = max(potential(*,loc,ialt,i)) - min(potential(*,loc,ialt,i))

cpcp = cpcp/1000.0

  plot, ut,cpcp, ytitle = 'kV', pos = pos, /noerase, $
	xtitle = 'UT Hours'

  oplot, [ut(it(0)),ut(it(0))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(1)),ut(it(1))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(2)),ut(it(2))], [0.0,10000.0], thick = 0.1

  xyouts, 0.1,175.0,' Cross Polar Cap Potential', charsize = 0.75

closedevice

mr = 60

loc = where(la lt mr)

setdevice, 'tn.ps', 'p', 4, 0.90

ppp = 6
space = 0.02
pos_space, ppp,space, sizes, ny = 2

readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

plotdumb

ialt = 5

get_position, ppp, space, sizes, 0, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pos_save = fltarr(4)
pos_save(0) = pos(0)
pos_save(2) = pos(2)

maxi = max(tn(*,loc,ialt,*))
mini = min(tn(*,loc,ialt,*))

levels = (maxi-mini)*findgen(9)/8.0 + mini
c_levels = (maxi-mini)*findgen(30)/29.0 + mini

it = [12,25,31]

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = tn(*,*,ialt,it(0))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /long

get_position, ppp, space, sizes, 1, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = tn(*,*,ialt,it(1))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long


xyouts, (pos(0)+pos(2))/2.0, pos(3) + 2*space, $
	'Neutral Temperature at '+tostr(fix(alts(ialt,nt/2)))+' Km', $
        charsize = 1.25,/norm,alignment = 0.5

get_position, ppp, space, sizes, 2, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = tn(*,*,ialt,it(2))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'K', /right, color = 0

get_position, ppp, space, sizes, 3, pos
pos_save(0) = pos(0) - space/2.0

get_position, ppp, space, sizes, 5, pos
pos_save(2) = pos(2) - space/2.0

pos_save(1) = pos(1)
pos_save(3) = pos(3)

pos = pos_save

cpcp = fltarr(nt)

for i=0,nt-1 do cpcp(i) = max(tn(*,loc,ialt,i))

  plot, ut,cpcp, ytitle = 'K', pos = pos, /noerase, $
	xtitle = 'UT Hours', yrange = mm(cpcp)

  oplot, [ut(it(0)),ut(it(0))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(1)),ut(it(1))], [0.0,10000.0], thick = 0.1
  oplot, [ut(it(2)),ut(it(2))], [0.0,10000.0], thick = 0.1

  xyouts, 0.1,403.0,'Maximum Temperature', charsize = 0.75

closedevice


setdevice, 'electron.ps', 'p', 4, 0.90

ppp = 6
space = 0.02
pos_space, ppp,space, sizes, ny = 2

readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

plotdumb

ialt = 5

get_position, ppp, space, sizes, 0, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pos_save = fltarr(4)
pos_save(0) = pos(0)
pos_save(2) = pos(2)

maxi = max(electron(*,loc,ialt,*))
mini = min(electron(*,loc,ialt,*))

levels = (maxi-mini)*findgen(9)/8.0 + mini
c_levels = (maxi-mini)*findgen(30)/29.0 + mini

it = [12,25,31]

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = electron(*,*,ialt,it(0))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /long

get_position, ppp, space, sizes, 1, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = electron(*,*,ialt,it(1))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long


xyouts, (pos(0)+pos(2))/2.0, pos(3) + 2*space, $
	'Electron Density at '+tostr(fix(alts(ialt,nt/2)))+' Km', $
        charsize = 1.25,/norm,alignment = 0.5

get_position, ppp, space, sizes, 2, pos
pos(0) = pos(0) - space/2.0
pos(2) = pos(2) - space/2.0

pot = fltarr(nlon+1,nlat)
pot(0:nlon-1,*) = electron(*,*,ialt,it(2))
pot(nlon,*) = pot(0,*)

contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase, $
	levels = c_levels, c_colors = c_colors,/cell_fill
contour, pot(*,loc), x(*,loc), y(*,loc), /follow,$
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], yrange=[-mr,mr],$
	/noerase
plotmlt, mr, /no06, /no00, /no18, /long

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
plotct, ncolors, ctpos, mm(levels)/1.0e5, 'x10!E5!N /cm!E3!N', $
	/right, color = 0

get_position, ppp, space, sizes, 3, pos
pos_save(0) = pos(0) - space/2.0

get_position, ppp, space, sizes, 5, pos
pos_save(2) = pos(2) - space/2.0

pos_save(1) = pos(1)
pos_save(3) = pos(3)

pos = pos_save

cpcp = fltarr(nt)

for i=0,nt-1 do cpcp(i) = max(electron(*,loc,ialt,i))/1.0e5

  plot, ut,cpcp, ytitle = 'x10!E5!N /cm!E3!N', pos = pos, /noerase, $
	xtitle = 'UT Hours', yrange = mm(cpcp)

  oplot, [ut(it(0)),ut(it(0))], [0.0,1.0e6], thick = 0.1
  oplot, [ut(it(1)),ut(it(1))], [0.0,1.0e6], thick = 0.1
  oplot, [ut(it(2)),ut(it(2))], [0.0,1.0e6], thick = 0.1

  xyouts, 0.1,2.4,'Maximum Electon Density', charsize = 0.75

closedevice

end