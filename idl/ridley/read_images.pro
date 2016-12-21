pro get_mlat, fin, fsize, glat, glon, mlat, mlon

  glatr = 0.25*float(round(glat/0.25))
  if glon lt 0.0 then glon = glon + 360.0
  glonr = 0.25*float(round(glon/0.25))
  dlat = long(4.0*(90.0-glatr))
  dlon = long(4.0*glonr)
  total = dlat*(360*4) + dlon
  total = long(16)*total
  if (total le fsize-16) then begin
    point_lun,fin,total
    readu,fin,glat,glon,mlat,mlon
    mlon = mlon-12.0
    if mlon lt 0.0 then mlon = mlon + 24.0
  endif else begin
    mlat = -999.0
    mlon = -999.0
  endelse

  return

end

pro read_images, filelist, nskip, uttimes, inttimes, image_save, min_lat

  flist = findfile(filelist)

  nfiles = n_elements(flist)

  fin = 2
  close,fin
  openu,fin,"/d/ridley/d.data/d.polar/d.software/d.savefiles/cgm_bin.dat"
  a = fstat(fin)
  fsize = a.size
  dx = 0.25

  uttimes = strarr(nfiles/nskip+1)
  inttimes = fltarr(nfiles/nskip+1)
  d1      = 0.0
  d2      = 0.0
  d3      = 0.0
  time    = ''

  isize=200
  image_save = fltarr(nfiles/nskip+1,isize,isize)

  ni = 0

  for n=0,nfiles-1,nskip do begin

    print, 'Reading file ',flist(n),' ; File ',tostr(n+1),' of ',tostr(nfiles)
    openr,1,flist(n)
    xsd = 0
    ysd = 0
    readf,1, xsd,ysd
    readf,1, time
    uttimes(ni) = strmid(time,0,13)
    if strlen(time) gt 14 then 					$
      inttimes(ni) = float(strmid(time,13,strlen(time)-13)) $
    else inttimes(ni) = 0.0

    a = fstat(1)
    sif = a.size
    si5 = long(xsd)*long(ysd)*long(4*5)

    imag = fltarr(xsd,ysd)
    lat = fltarr(xsd,ysd)
    lon = fltarr(xsd,ysd)

    readu,1,lat
    readu,1,lon
    readu,1,imag

    cgm = 0
    if sif gt si5 then begin
      mlat = fltarr(xsd,ysd)
      mlt  = fltarr(xsd,ysd)
      readu,1,mlat
      readu,1,mlt
      cgm = 1
    endif else print, '  converting coordinates to cgm'

    close,1

    si = n_elements(image_save(0,0,*))
    fimage = fltarr(si,si)
    nimage = fltarr(si,si)

    ran = 90.0 - min_lat
    fac = float(si)/float(ran*2)
    offs = si/2
    ut = float(strmid(uttimes(ni),7,2))+float(strmid(uttimes(ni),9,2))/60.0

    for i=0,xsd-1 do for j=0,ysd-1 do begin

      if (imag(i,j) gt 0.0) then begin
        if (not cgm) then 					$
	  get_mlat, fin, fsize, lat(i,j),lon(i,j), ml, mm	$
	else begin
	  ml = mlat(i,j)
	  mm = mlt(i,j)
	endelse

        if (ml gt -24.0) and (ml gt min_lat) then begin

          clat = 90.0-ml
          the = !pi*(ut-mm)/12.0 + !pi/2.0
	  x = round(clat*cos(the)*fac)+offs
	  y = round(clat*sin(the)*fac)+offs
	  if (x gt 1) and			$
	     (x le si-2) and			$
	     (y gt 1) and			$
	     (y le si-2) then begin
	    fimage(x,y) = fimage(x,y) + imag(i,j)
	    nimage(x,y) = nimage(x,y) + 1.0
	    fimage(x-1,y) = fimage(x-1,y) + 0.25*imag(i,j)
	    nimage(x-1,y) = nimage(x-1,y) + 0.25
	    fimage(x,y-1) = fimage(x,y-1) + 0.25*imag(i,j)
	    nimage(x,y-1) = nimage(x,y-1) + 0.25
	    fimage(x+1,y) = fimage(x+1,y) + 0.25*imag(i,j)
	    nimage(x+1,y) = nimage(x+1,y) + 0.25
	    fimage(x,y+1) = fimage(x,y+1) + 0.25*imag(i,j)
	    nimage(x,y+1) = nimage(x,y+1) + 0.25
          endif

        endif

      endif

    endfor

    loc = where(nimage ne 0.0,count)
    if count gt 0 then fimage(loc) = fimage(loc)/float(nimage(loc))

    image_save(ni,*,*) = fimage

    ni = ni + 1

  endfor

  close, fin

  uttimes = uttimes(0:ni-1)
  image_save = image_save(0:ni-1,*,*)

  return

end
