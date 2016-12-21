
setdevice

dlat = 0.25
dlon = 1.0
minlat = 50.0
missing = -1.0e32

flist = ''
read, 'Enter file list to read : ', flist

nskip = 0
read, 'enter number to skip between each image (0 for none) : ', nskip
nskip = nskip + 1

read_raw_image,flist,nskip,uttimes, inttimes, image_save, coor_save

; geolt = (ut + geolon/15.0) mod 24.0

npts = n_elements(uttimes)

i = long(0)

nlat = (90.0-minlat)/dlat+1
nlon = 360.0/dlon+1

boundary = fltarr(npts,2,nlon)

ppp = npts
pos_space,ppp,0.01,sizes

for i=long(0),npts-1 do begin

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
    aves = (lin_image(j,0:nlat-1))^2.0
    locnz = where(aves gt 0.0, countnz)
    if countnz gt 0 then begin
      mn = mean(aves(locnz))
      loc2 = where(aves gt mn,countk)
      if countk gt 0 then begin
        blat_p(j) = 90.0 - float(loc2(0))*dlat
        blat_e(j) = 90.0 - float(loc2(countk-1))*dlat
      endif
    endif
  endfor

  blat_p = [blat_p(nlon-10:nlon-1),blat_p,blat_p(0:9)]
  run_ave,blat_p,20,0.0
  blat_p = blat_p(10:n_elements(blat_p)-11)

  blat_e = [blat_e(nlon-10:nlon-1),blat_e,blat_e(0:9)]
  run_ave,blat_e,20,0.0
  blat_e = blat_e(10:n_elements(blat_e)-11)

  boundary(i,0,*) = blat_p
  boundary(i,1,*) = blat_e

  get_position,ppp,0.01,sizes,i,pos

  image1 = 255.0*image1/500.0
  image1 = 255.0 - image1

  tv, image1, pos(0), pos(1),	 				$
    xsize=(pos(2)-pos(0)), ysize=(pos(3)-pos(1)),		$
    /norm

  plot, blon, 90.0-blat_p,/noerase, 				$
	xrange = [0,360], yrange = [0,90.0-minlat],		$
	pos = pos, xstyle =1, ystyle = 1
  oplot, blon, 90.0-blat_e

endfor

plotdumb

for i=long(0),npts-1 do begin

  ut = float(strmid(uttimes(i),7,2)) + float(strmid(uttimes(i),9,2))/60.0

  blat_p(0:nlon-1) = boundary(i,0,*)
  blat_e(0:nlon-1) = boundary(i,1,*) 

  geolt = (ut + blon/15.0) mod 24.0
  the = (geolt-12.0)*!pi/12.0 + !pi/2.0
  ran = 90.0-blat_p

  get_position,ppp,0.01,sizes,i,pos

  plot, ran*cos(the),ran*sin(the), xrange=[-40,40],yrange=[-40,40],	$
	/noerase, pos = pos, xstyle = 5, ystyle = 5

  ran = 90.0-blat_e
  oplot, ran*cos(the),ran*sin(the)

  plotmlt, 40.0, /no00, /no06, /no12, /no18

endfor

if !d.name eq 'PS' then begin
  device, /close
  set_plot, 'X'
endif

end
