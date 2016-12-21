
loadct, 5

dlat = 0.5
dlon = 2.0
minlat = 50.0
missing = -1.0e32

flist = ''
read, 'Enter file list to read : ', flist

nskip = 0
;read, 'enter number to skip between each image (0 for none) : ', nskip
nskip = nskip + 1

read_raw_image,flist,nskip,uttimes, inttimes, image_save, coor_save

npts = n_elements(uttimes)

i = long(0)

nlat = (90.0-minlat)/dlat+1
nlon = 360.0/dlon+1

ppp = 6
pos_space,ppp,0.05,sizes

area = fltarr(npts)
uts = fltarr(npts)
re = 6270.0

for i=long(0),npts-1 do begin

  print, 'Working on image ',tostr(fix(i))

  pn = i mod ppp

  ut = float(strmid(uttimes(i),7,2)) + float(strmid(uttimes(i),9,2))/60.0
  rotation = (12.0-ut)*180.0/12.0 + 180.0
  if rotation lt 0.0 then rotation = rotation + 360.0

  lin_image = fltarr(nlon,nlat);+missing

  for j=0,n_elements(image_save(i,*,0))-1 do 				$
    for k=0,n_elements(image_save(i,0,*))-1 do begin
;      coor_save(i,1,j,k) = coor_save(i,1,j,k) + rotation
;      if coor_save(i,1,j,k) gt 360.0 then 				$
;	coor_save(i,1,j,k) = coor_save(i,1,j,k) - 360.0
      x = coor_save(i,1,j,k)/dlon
      y = (90.0-coor_save(i,0,j,k))/dlat
      if x ge 0 and x lt nlon and y ge 0 and y lt nlat then 		$
        lin_image(x,y) = image_save(i,j,k)
    endfor

;  smooth_image, lin_image, missing
;  expand_array, lin_image, 3

  nlat = n_elements(lin_image(0,*))
  nlon = n_elements(lin_image(*,0))

  blat_e = fltarr(nlon)
  blat_p = fltarr(nlon)
  blon = findgen(nlon)*dlon

  image1 = lin_image

  locnz = where(lin_image gt 0.0, countnz)
  if countnz gt 0 then 							$
    lin_image(locnz) = lin_image(locnz) - mean(lin_image(locnz))/2.0
  locnz = where(lin_image lt 0.0, countnz)
  if countnz gt 0 then lin_image(locnz) = 0.0
  lin_image = lin_image^2.0

  for j=0,n_elements(lin_image(*,0))-1 do begin
    loc = where(lin_image(j,*) eq max(lin_image(j,*)))
    aves = (lin_image(j,0:nlat-1)) ; ^2.0
    locnz = where(aves gt 0.0, countnz)
    if countnz gt 0 then begin
      mn = 0.75*mean(aves(locnz))
;      loc2 = where(aves gt mn,countk)
;      if countk gt 0 then begin
      locp = -1
      n = loc(0)
      while (locp eq -1) do 						$
	if (aves(n) lt mn) or (n eq 0) then locp = n else n=n-1
      loce = -1
      n = loc(0)
      while (loce eq -1) do 						$
	if (aves(n) lt mn) or (n eq nlat-1) then loce = n else n=n+1

        blat_p(j) = 90.0 - float(locp)*dlat
        blat_e(j) = 90.0 - float(loce)*dlat
;      endif
    endif
  endfor

  blat_p = [blat_p(nlon-10:nlon-1),blat_p,blat_p(0:9)]
  run_ave,blat_p,10,0.0
  blat_p = blat_p(10:n_elements(blat_p)-11)

  blat_e = [blat_e(nlon-10:nlon-1),blat_e,blat_e(0:9)]
  run_ave,blat_e,10,0.0
  blat_e = blat_e(10:n_elements(blat_e)-11)


  image2 = lin_image



  si = 250
  fac = float(si)/(2.0*(90.0-minlat))

  fimage = fltarr(si,si)
  nimage = fltarr(si,si)

;  for j=0,n_elements(image_save(i,*,0))-1 do 				$
;    for k=0,n_elements(image_save(i,0,*))-1 do begin
;  for j=0,n_elements(lin_image(0,*))-1 do 				$
;    for k=0,n_elements(lin_image(*,0))-1 do begin

;      kk = fix(coor_save(i,1,j,k))
;      jj = fix((90.0-coor_save(i,0,j,k))/dlat)

;      clat = fac*(90.0-coor_save(i,0,j,k))
;          the = !pi*(ut-mm)/12.0 + !pi/2.0
;      the = !pi*coor_save(i,1,j,k)/180.0

;      clat = fac*(float(j)*dlat)
;      the = !pi*float(k)*dlon/180.0

;      x = round(clat*cos(the)+float(si)/2.0)
;      y = round(clat*sin(the)+float(si)/2.0)
;      if (x gt 1) and			$
;	 (x le si-2) and			$
;	 (y gt 1) and			$
;	 (y le si-2) then begin
;	fimage(x,y) = fimage(x,y) + lin_image(k,j)
;	nimage(x,y) = nimage(x,y) + 1.0
;	fimage(x-1,y) = fimage(x-1,y) + 0.25*lin_image(k,j)
;	nimage(x-1,y) = nimage(x-1,y) + 0.25
;	fimage(x,y-1) = fimage(x,y-1) + 0.25*lin_image(k,j)
;	nimage(x,y-1) = nimage(x,y-1) + 0.25
;	fimage(x+1,y) = fimage(x+1,y) + 0.25*lin_image(k,j)
;	nimage(x+1,y) = nimage(x+1,y) + 0.25
;	fimage(x,y+1) = fimage(x,y+1) + 0.25*lin_image(k,j)
;	nimage(x,y+1) = nimage(x,y+1) + 0.25
;	fimage(x,y) = fimage(x,y) + image_save(i,j,k)
;	nimage(x,y) = nimage(x,y) + 1.0
;	fimage(x-1,y) = fimage(x-1,y) + 0.25*image_save(i,j,k)
;	nimage(x-1,y) = nimage(x-1,y) + 0.25
;	fimage(x,y-1) = fimage(x,y-1) + 0.25*image_save(i,j,k)
;	nimage(x,y-1) = nimage(x,y-1) + 0.25
;	fimage(x+1,y) = fimage(x+1,y) + 0.25*image_save(i,j,k)
;	nimage(x+1,y) = nimage(x+1,y) + 0.25
;	fimage(x,y+1) = fimage(x,y+1) + 0.25*image_save(i,j,k)
;	nimage(x,y+1) = nimage(x,y+1) + 0.25
;      endif

;    endfor


  loc = where(nimage gt 0,count)
  if count gt 0 then fimage(loc) = fimage(loc)/nimage(loc)
;  loc = where(nimage eq 0,count)
;  if count gt 0 then fimage(loc) = missing

;  smooth_image, fimage, missing
;  loc = where(fimage eq missing,count)
;  if count gt 0 then fimage(loc) = 0.0

;  flat = [90.0,80.0,70.0,60.0]
;  flon = [0.0,00.0,00.0,00.0]
;  clat = fac*(90.0-flat)
;  the = !pi*flon/180.0
;  x = round(clat*cos(the)+float(si)/2.0)
;  y = round(clat*sin(the)+float(si)/2.0)
;  fimage(x,y) = max(image_save)

  if pn eq 0 then plotdumb

  get_position,ppp,0.05,sizes,pn,pos

  image1 = 255.0*image1/1000.0
;  lin_image = 255.0*lin_image/400.0

;  image1 = max(image1) - image1

  tv, image1, pos(0), pos(1), 				$
    xsize=(pos(2)-pos(0)), ysize=(pos(3)-pos(1)),		$
    /norm

  plot, blon, 90.0-blat_p,/noerase, 				$
	xrange = [0,360], yrange = [0,90.0-minlat],		$
	pos = pos, xstyle =1, ystyle = 1
  oplot, blon, 90.0-blat_p, color=255
  oplot, blon, 90.0-blat_e, color=255

  area(i) = 0

  dl = (blon(1) - blon(0))*!pi/180.0

  for j=0,n_elements(blat_p)-1 do 				$
    area(i) = area(i) + dl*(re*((90.0-blat_p(j))*!pi/180.0))^2.0

  uts(i) = float(strmid(uttimes(i),7,2)) + 			$
	   float(strmid(uttimes(i),9,2))/60.0 +			$
	   float(strmid(uttimes(i),9,2))/3600.0

;  get_position,1,0.05,sizes,1,pos

;  fimage = max(fimage) - fimage

;  lat = 90.0
;  long = 0.0
;  minlat = 50
;  maxlat = 90
;  minlon = 0
;  maxlon = 360
;  limit = [minlat, minlon, maxlat, maxlon]
;  !p.region = pos
;  map_set,lat, long, rotation, limit = limit, 			$
;	/stereo, /noborder, /advance, /cont

;  i_lat = fltarr(5)
;  i_lon = fltarr(5)

;  maxval = max(image_save)

;  nx = n_elements(image_save(i,*,0))
;  ny = n_elements(image_save(i,0,*))

;  for j=1,nx-2 do for k=1,ny-2 do begin

;    np = 0
;    tp = 0
;    for l=j-1,j+1,2 do for m=k-1,k+1,2 do			$
;      if (coor_save(i,1,l,m) ge 360+minlon) then			$
;	coor_save(i,1,l,m) = coor_save(i,1,l,m) - 360.0
;    for l=j-1,j+1,2 do for m=k-1,k+1,2 do begin
;      if (coor_save(i,0,l,m) ge minlat) and			$
;	 (coor_save(i,0,l,m) le maxlat) and			$
;	 (coor_save(i,1,l,m) ge minlon) and			$
;	 (coor_save(i,1,l,m) le maxlon) then np = np + 1
;      tp = tp+1
;    endfor

;    if tp eq np then begin

;	i_lat(0) = coor_save(i,0,j-1,k-1)
;	i_lat(1) = coor_save(i,0,j-1,k+1)
;	i_lat(2) = coor_save(i,0,j+1,k+1)
;	i_lat(3) = coor_save(i,0,j+1,k-1)
;	i_lat(4) = coor_save(i,0,j-1,k-1)
;	i_lon(0) = coor_save(i,1,j-1,k-1)
;	i_lon(1) = coor_save(i,1,j-1,k+1)
;	i_lon(2) = coor_save(i,1,j+1,k+1)
;	i_lon(3) = coor_save(i,1,j+1,k-1)
;	i_lon(4) = coor_save(i,1,j-1,k-1)
;        colo = 255.0*image_save(i,j,k)/maxval

;        polyfill, i_lon, i_lat, color=colo
;    endif
;  endfor

;    MAP_CONTINENTS,MLINETHICK=1, color = 255

;  tvscl, fimage, pos(0), pos(1), 				$
;    xsize=(pos(2)-pos(0)), ysize=(pos(3)-pos(1)),		$
;    /norm

;  plot, (90.0-blat_p)*cos(blon*!pi/180.0),			$
;	(90.0-blat_p)*sin(blon*!pi/180.0),/noerase, 		$
;	xrange = [-(90.0-minlat),(90.0-minlat)], 		$
;	yrange = [-(90.0-minlat),(90.0-minlat)],		$
;	xstyle = 1, ystyle = 1, pos = pos

  oplot, blon, blat_p, color = 255
  oplot, blon, blat_e, color = 255

;  oplot, (90.0-blat_e)*cos(blon*!pi/180.0),			$
;	 (90.0-blat_e)*sin(blon*!pi/180.0)

;  oplot, [0.0,0.0],[-40.0,40.0]
;  oplot, [10.0,10.0],[-40.0,40.0]
;  oplot, [20.0,20.0],[-40.0,40.0]
;  oplot, [30.0,30.0],[-40.0,40.0]
;  oplot, [-40.0,40.0], [0.0,0.0]

;  if i ne npts-1 then prompt_for_next

endfor

run_ave,area,5,0.0

plot,uts,area, ystyle = 1, xstyle = 1

if !d.name eq 'PS' then begin
  device, /close
  set_plot, 'X'
endif

end
