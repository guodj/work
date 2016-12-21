minlat = 40.0
nlats  = (90.0-minlat)/5.0
nlongs = 24

filein = ask('data file from T-GCM run','NO_eq.bin')
bora = ''
read,'Is this file a binary or ascii file (b/a) [b] : ',bora

psfile = ask('PS file name (return for screen)','')
if strlen(psfile) eq 0 then window else setdevice, psfile, 'p'

print, ''
print, '1. Geographic Latitude vs. Geographic Longitude (1 height)'
print, '2. Magnetic Latitude vs. Magnetic Local Time (1 height)'
plottype = fix(ask('type of plot','1'))

if (strmid(bora,0,1) eq 'b') or (strmid(bora,0,1) eq 'B') then 		$
  rd_tgcm_bin, filein, in_data, vars, utdes, axdes, xax, yax		$
else									$
  rd_tgcm_acsii, filein, in_data, vars, utdes, axdes, xax, yax, 200, run

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

for iut = 0, nut-2 do begin

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

if strpos(vars(0),'DENSITY') eq -1 then data = 10.0^data

loc = where(data ne 1.0,count)

if count gt 0 then begin

  data = data(loc)
  latdata = latdata(loc)
  londata = londata(loc)
  mltdata = mltdata(loc)
  mladata = mladata(loc)

endif

mini = 1.05*min(data)
maxi = 0.90*max(data)

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

outbin = fltarr(nlongs,nlats)

maxrange = 90.0 - minlat

for j=0,1 do begin

  if j eq 0 then begin
    plotdumb
    xyouts, 0.0, 0.0, 'TIMEGCM Equinox Condition', /norm
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
      if count gt 0 then outbin(ix,iy) = mean(data(loc))
    endfor
    print, ix+1, ' out of ',nlongs
  endfor

  y = outbin*0.0 + 1.0
  smooth, outbin, y, 0.0

  plot, [-maxrange,maxrange],[-maxrange,maxrange],		$
	xstyle = 5, ystyle = 5, pos = pos,			$
	/nodata, /noerase

  xcors = outbin*0.0
  ycors = outbin*0.0
  datai = outbin

  image = ncolors*(outbin-mini)/(maxi-mini)

  for ix=0,nlongs-1 do begin
    lon1 = float(ix)*dlong
    for iy=0,nlats-1 do begin
      lat1 = maxrange - float(iy)*dlat
      lat2 = maxrange - float(iy+1)*dlat
      lats = [fltarr(11)+lat1, fltarr(11)+lat2, lat1]
      lons = [findgen(11)*dlong/10.0 + lon1,			$
	      (10.0-findgen(11))*dlong/10.0 + lon1,		$
              lon1] - 90.0
      x = lats*cos(lons*!pi/180.0)
      y = lats*sin(lons*!pi/180.0)
      polyfill, x, y, color = image(ix,iy)
      xcors(ix,iy) = mean(x)
      ycors(ix,iy) = mean(y)
    endfor
  endfor

  contour, datai, xcors, ycors, /follow, nlevels = 20, 		$
           color = 0, /noerase, pos = pos, xstyle = 5, 		$
           ystyle = 5, xrange = [-maxrange,maxrange],		$
           yrange = [-maxrange,maxrange]

  if plottype eq 1 then 					$
    plotmlt, maxrange, /longs					$
  else								$
    plotmlt, maxrange

  plotmlt, maxrange, /no00, /no06, /no12, /no18, /white
  xyouts, pos(0), pos(3)+0.01, title, /norm

  posct = [0.90,pos(1),0.92,pos(3)]
  plotct, ncolors, posct, [mini,maxi], strcompress(vars(0)), /right

endfor

closedevice

end