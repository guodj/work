minlat = 40.0
nlats  = (90.0-minlat)/5.0
nlongs = 48

filein = ask('data file from T-GCM run','NO_eq.bin')
bora = ''
read,'Is this file a binary or ascii file (b/a) [a] : ',bora

psfile = ask('PS file name (return for screen)','')
if strlen(psfile) eq 0 then window else setdevice, psfile, 'p'

print, ''
print, '1. Geographic Latitude vs. Geographic Longitude (1 height)'
print, '2. Magnetic Latitude vs. Magnetic Local Time (1 height)'
plottype = fix(ask('type of plot','1'))

timegcm = ask('title','TIMEGCM Output')

if (strmid(bora,0,1) eq 'b') or (strmid(bora,0,1) eq 'B') then 		$
  rd_tgcm_bin, filein, in_data, vars, utdes, axdes, xax, yax		$
else									$
  rd_tgcm_ascii, filein, in_data, vars, utdes, axdes, xax, yax, 200, run

print, '1. Plot normal data'
print, '2. Plot log(data)'
print, '3. Plot 10^data'
plotdata = ask('Plotting option (1,2,3)','1')

if (plotdata eq '2') then in_data = alog10(in_data)
if (plotdata eq '3') then in_data = 10^(in_data)

xax = (xax+360) mod 360.0

npts = n_elements(in_data)

data    = fltarr(npts)
latdata = fltarr(npts)
londata = fltarr(npts)
mladata = fltarr(npts)
mlodata = fltarr(npts)
mltdata = fltarr(npts)
time    = dblarr(npts)

itime = [1995,03,21,0,0,0]
c_a_to_r, itime, basetime
utadd = float(strmid(utdes,3,4))*3600.0

itime = [1995,01,01,0,0,0]
c_a_to_r, itime, bt_year

nut = n_elements(in_data(*,0,0))
nlo = n_elements(in_data(0,*,0))
nla = n_elements(in_data(0,0,*))

mlat = fltarr(nla,nlo)
mlon = fltarr(nla,nlo)

year = itime(0)

for ilo = 0, nlo-1 do begin

  for ila = 0, nla-1 do begin

    lon = xax(ilo)
    lat = yax(ila)
    mag_pos   = cnvcoord(lat,lon,120.0)
    mlat(ila,ilo) = mag_pos(0)
    mlon(ila,ilo) = mag_pos(1)

  endfor

endfor

for iut = 0, nut-1 do begin

  for ilo = 0, nlo-1 do begin

    for ila = 0, nla-1 do begin

      index = long(ila) + long(ilo*nla) + long(iut)*long(nla*nlo)
      data(index) = in_data(iut,ilo,ila)
      londata(index) = xax(ilo)
      latdata(index) = yax(ila)
      mladata(index) = mlat(ila,ilo)
      mlodata(index) = mlon(ila,ilo)
      if plottype eq 2 then begin
        time(index)    = basetime + utadd(iut)
        utsec          = time(index) - bt_year
        mag_mlt        = mlt(year, long(utsec), mlodata(index))
        mltdata(index) = mag_mlt(0)
      endif

    endfor

  endfor

  print, iut+1,' out of ',nut-1

endfor

loc = where(data ne 1.0,count)

if count gt 0 then begin

  data = data(loc)
  latdata = latdata(loc)
  londata = londata(loc)
  mltdata = mltdata(loc)
  mladata = mladata(loc)

endif

mini = min(data)
maxi = max(data)

if plottype eq 2 then begin

  latdata = mladata
  londata = mltdata*360.0/24.0
  h_string = 'Height: 120 km ; Magnetic Coordinates'

endif else begin

  h_string = 'Height: 120 km ; Geographic Coordinates'

endelse

; read in color table

ct_dir = getenv('IDL_EXTRAS')
ctname = ct_dir+'joule4.ct'
readct,ncolors,ctname

dx = 1.25*float(!d.y_ch_size)/float(!d.y_size)
ppp = 2
space = dx*2.0
pos_space,ppp,space,sizes, ny = ppp

outbin = fltarr(nlongs+1,nlats+1)
xcors = outbin*0.0
ycors = outbin*0.0
lacor = outbin*0.0
locor = outbin*0.0

maxrange = 90.0 - minlat

for j=0,1 do begin

  if j eq 0 then begin
    plotdumb
    xyouts, 0.0, 0.0, timegcm, /norm
    xyouts, 0.0, -dx*1.5, h_string, /norm
  endif

  if (j eq 0) then title = 'Northern Hemisphere'
  if (j eq 1) then title = 'Southern Hemisphere'

  get_position, ppp, space, sizes, j mod ppp, pos

  dlong = 360.0/(nlongs)
  dlat  = (90.0-minlat)/(nlats)
  if j eq 0 then fac = 1.0 else fac = -1.0

  outbin(*,*) = 0.0

  for ix=0,nlongs-1 do begin
    lon1 = float(ix)*dlong
    lon2 = float(ix+1)*dlong
    for iy=0,nlats-1 do begin
      lat1 = minlat + float(iy)*dlat
      lat2 = minlat + float(iy+1)*dlat
      loc = where( (fac*latdata ge lat1) and			$
	           (fac*latdata lt lat2) and			$
		   (londata ge lon1) and			$
		   (londata lt lon2),count)
      if count gt 0 then outbin(ix,iy) = amean(data(loc))
      x = (90.0-(lat1+lat2)/2.0)*cos(((lon1+lon2)/2.0)*!pi/180.0)
      y = (90.0-(lat1+lat2)/2.0)*sin(((lon1+lon2)/2.0)*!pi/180.0)
      xcors(ix,iy) = x
      ycors(ix,iy) = y
      lacor(ix,iy) = (lat1+lat2)/2.0
      locor(ix,iy) = (lon1+lon2)/2.0
    endfor
  endfor

  outbin(nlongs,*) = outbin(0,*)
  outbin(*,nlats) = mean(outbin(*,nlats-1))

  y = outbin*0.0 + 1.0
;  smooth, outbin, y, 0.0

  plot, [-maxrange,maxrange],[-maxrange,maxrange],		$
	xstyle = 5, ystyle = 5, pos = pos,			$
	/nodata, /noerase

  datai = outbin

  image = ncolors*(outbin-mini)/(maxi-mini)

  xcors(nlongs,*) = xcors(0,*)
  ycors(nlongs,*) = ycors(0,*)
  xcors(*,nlats) = 0.0
  ycors(*,nlats) = 0.0

  lacor(nlongs,*) = lacor(0,*)
  locor(nlongs,*) = locor(0,*)
  lacor(*,nlats) = 90.0
  locor(*,nlats) = locor(*,0)

  levels = float(ncolors)*findgen(30)/29

  contour, image, xcors, ycors, nlevels=30, 		$
           /noerase, pos = pos, xstyle = 5, 		$
           ystyle = 5, xrange = [-maxrange,maxrange],		$
           yrange = [-maxrange,maxrange],/fill,levels=levels, c_colors=levels

  contour, datai, xcors, ycors, /follow, nlevels=15, 		$
           /noerase, pos = pos, xstyle = 5, 		$
           ystyle = 5, xrange = [-maxrange,maxrange],		$
           yrange = [-maxrange,maxrange]

  if plottype eq 1 then 					$
    plotmlt, maxrange, /longs					$
  else								$
    plotmlt, maxrange

  xyouts, pos(0), pos(3)+0.01, title, /norm

  posct = [0.90,pos(1),0.92,pos(3)]
  plotct, ncolors, posct, [mini,maxi], strcompress(vars(0)), /right,color=0

endfor

closedevice

end