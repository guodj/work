pro read_block, flatfile, times, data

  Setup_Var_Menu, 1, flatfile, timeint, vars, units, ncol, nrows, rowlen
  nr = nrows(0)
  nc = ncol(0)+1

  data = fltarr(nc,nr)
  time = dblarr(nr)
  dumb = fltarr(nc)
  dtime = double(0.0)

  openr, 11, flatfile+'.dat'

  i = long(0)

  for i=long(0),long(nr-1) do begin

    readu, 11, dtime
    readu, 11, dumb

    time(i) = dtime
    data(*,i) = dumb

  endfor

  close,11

  return

end

pro plotorbit, filelist, psfile, nbin, nskip

  flist = findfile(filelist)
  nfile = n_elements(flist)

  trow = long(0)

  for i=0,nfile-1 do begin

    datafile = flist(i)
    loc = strpos(datafile, ".hed")
    datafile = strmid(datafile,0,loc)
    flist(i) = datafile
    Setup_Var_Menu, 1, datafile, timeint, vars, units, ncol, nrows, rowlen
    trow = trow + nrows(0)

  endfor

  print, 'total rows : ',trow
  print, 'ncol :',ncol(0)

  data = fltarr(ncol(0)+1,trow)
  height = fltarr(nfile,2)

  crow = 0

  for i=0,nfile-1 do begin

    datafile = flist(i)
    print, 'reading file : ',datafile

    read_block, datafile, time, data_temp

    ctrow = n_elements(data_temp(0,*))
    data(*,crow:crow+ctrow-1) = data_temp(*,0:ctrow-1)

    height(i,*) = mm(data_temp(11,*))
    crow = crow + ctrow

  endfor

  minlat = 50.0

  maxran = 90.0-minlat

  lat = data(17,*)
  mlt = data(15,*)
  alt = data(11,*)

  loc = where(lat gt minlat,count)

  if strlen(psfile) gt 0 then setdevice,psfile else plotdumb

  ppp = nbin+1
  space = 0.05
  pos_space,ppp,space,sizes

  dz = 1000.0/float(nbin)

  get_position,ppp,space,sizes,0,pos

  yr = mm(height)
  plot, height(*,0), xstyle = 1, yrange = yr, pos = pos,		$
	ytitle = 'Altitude (km)'
  oplot, height(*,1)

  if count gt 0 then begin

    for i=0,nbin-1 do begin

      get_position,ppp,space,sizes,i+1,pos

      plot, [-maxran,maxran],[-maxran,maxran], xstyle=5,ystyle=5,	$
	pos=pos,/nodata, /noerase
      xyouts,-maxran,-maxran, tostr((i+1)*dz)+' km',charsize=0.75
      plotmlt,maxran

      l2 = where(alt(loc) gt float(i)*dz and 				$
		   alt(loc) le float(i+1)*dz, count)

      xyouts,maxran,-maxran, tostrf(float(count)),charsize=0.75, 	$
	alignment=1.0

      if count gt 0 then begin

	j=long(0)

	for j=long(0),long(count-1),nskip do begin

	  ind = loc(l2(j))

          x = (90.0-lat(ind))*cos(mlt(ind)*!pi/12.0 - !pi/2.0)
          y = (90.0-lat(ind))*sin(mlt(ind)*!pi/12.0 - !pi/2.0)

          oplot, [x],[y], psym = 3

	endfor

      endif

    endfor

  endif

  if !d.name eq 'PS' then begin
    device, /close
    set_plot,'X'
  endif

  return

end

print, 'Enter file list to print : '
flist = ''
read,flist

print, 'Enter ps file name (return for screen) : '
psfile = ''
read,psfile

print, 'Number of altitude bins between 0 - 1000 km : '
nbin = 0
read,nbin
if nbin le 0 then nbin = 1

print, 'Number of pixels to skip between each plotted (99) : '
nskip = 0
read,nskip
nskip = nskip+1

plotorbit, flist, psfile, nbin, nskip

end
