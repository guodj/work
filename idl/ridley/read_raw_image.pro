pro read_raw_image, filelist, nskip, uttimes, inttimes, image_save, coor_save

  flist = findfile(filelist)

  nfiles = n_elements(flist)

  uttimes = strarr(nfiles/nskip+1)
  inttimes = fltarr(nfiles/nskip+1)
  d1      = 0.0
  d2      = 0.0
  d3      = 0.0
  time    = ''

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

    if n eq 0 then begin
      image_save = fltarr(nfiles/nskip+1,xsd,ysd)
      coor_save = fltarr(nfiles/nskip+1,4,xsd,ysd)
    endif

    imag = fltarr(xsd,ysd)
    lat = fltarr(xsd,ysd)
    lon = fltarr(xsd,ysd)

    readu,1,lat
    readu,1,lon
    readu,1,imag

    image_save(ni,*,*) = imag(*,*)
    coor_save(ni,0,*,*) = lat(*,*)
    coor_save(ni,1,*,*) = lon(*,*)

    if sif gt si5 then begin
      mlat = fltarr(xsd,ysd)
      mlt  = fltarr(xsd,ysd)
      readu,1,mlat
      readu,1,mlt
      coor_save(ni,2,*,*) = mlat(*,*)
      coor_save(ni,3,*,*) = mlt(*,*)
    endif else begin
      coor_save(ni,2,*,*) = -1.0e32
      coor_save(ni,3,*,*) = -1.0e32
    endelse

    close,1

    ni = ni + 1

  endfor

  uttimes = uttimes(0:ni-1)
  image_save = image_save(0:ni-1,*,*)
  coor_save = coor_save(0:ni-1,0:3,*,*)

  return

end
