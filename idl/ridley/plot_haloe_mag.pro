
common ffinfo, header

; ----------------------------------------------------------------
; BIG ASSUMPTIONS :

local_time_id  = [0]
local_time_col = [2]

lat_col = [0,1]
lat_id  = intarr(n_elements(lat_col))*0

heights = [1.102e-01, 9.127e-01, 6.125e+00, 8.822e+00, 1.133e+01, 1.378e+01, $
           1.622e+01, 1.867e+01, 2.112e+01, 2.360e+01, 2.610e+01, 2.864e+01, $
           3.121e+01, 3.382e+01, 3.651e+01, 3.929e+01, 4.216e+01, 4.512e+01, $
           4.818e+01, 5.127e+01, 5.432e+01, 5.727e+01, 6.014e+01, 6.292e+01, $
           6.562e+01, 6.823e+01, 7.077e+01, 7.323e+01, 7.564e+01, 7.800e+01, $
           8.030e+01, 8.256e+01, 8.477e+01, 8.693e+01, 8.910e+01, 9.126e+01, $
           9.344e+01, 9.562e+01, 9.783e+01, 1.001e+02, 1.024e+02, 1.049e+02, $
           1.075e+02, 1.103e+02, 1.136e+02, 1.175e+02, 1.221e+02, 1.275e+02, $
           1.338e+02, 1.411e+02, 1.495e+02, 1.590e+02, 1.695e+02, 1.813e+02, $
           1.941e+02]

; ----------------------------------------------------------------

haloefile = ''
read, 'Enter Flatfile name for HALOE data : ', haloefile

paramfile = ''
read, 'Enter Flatfile name for index data (return for none) : ', paramfile

psfile = ''
read, 'Enter ps filename (return for screen) : ',psfile

if strlen(psfile) gt 0 then setdevice, psfile,'p',4,0.90		$
else window,0

nfile = 1
datafile = strarr(6)
datafile(0) = haloefile
datafile(1) = paramfile
if strlen(paramfile) gt 0 then nfile=2

; get variable names:

Setup_Var_Menu, nfile, datafile, timeint, titles, units,	$
                count, nrows, rowlen, colloc, coltype, colnele

header = {nf : nfile, df : datafile, ti : timeint,		$  
   	  na : titles, un : units, nr : nrows, rl : rowlen,	$
          loc : colloc, type : coltype, nele : colnele,		$
   	  unit : intarr(nfile+1), ncol : count+1}

; ask which variables to plot:

print, ''
for id=0,nfile-1 do for col = 0,header.ncol(id)-1 do 			$
  print, tostr(id+1),', ',chopr('000'+tostr(col+1),2),'. ',header.na(id,col)

done = 0
nvar = 0

while not done do begin

  cid = 0
  ccol = 0
  read, 'Enter file, column to plot (0,0 to end entering) : ',cid,ccol

  if (cid gt 0) and (ccol gt 0) then begin

    if (nvar eq 0) then begin

      id_list = cid-1
      col_list = ccol-1

    endif else begin

      id_list  = [id_list, cid-1]
      col_list = [col_list,ccol-1]

    endelse

    nvar = nvar + 1

  endif else done = 1

endwhile

print, '1. No Contours on color plots'
print, '2. Contours on color plots'

con = fix(ask('selection','1'))

if nvar gt 0 then begin

; break apart data into scalor and vector quantities:

  id  = intarr(nvar,2)
  col = intarr(nvar,2)

  nvar_scal = 0
  nvar_vect = 0

  for i=0,nvar-1 do 						$
    if (header.nele(id_list(i),col_list(i),0) eq 1) then begin
      id(nvar_scal,0) = id_list(i)
      col(nvar_scal,0) = col_list(i)
      nvar_scal = nvar_scal + 1
    endif else begin
      id(nvar_vect,1) = id_list(i)
      col(nvar_vect,1) = col_list(i)
      nvar_vect = nvar_vect + 1
    endelse

  print, ''

  sstime = ask('Starting Time to plot',header.ti(0,0))
  setime = ask('Ending Time to plot  ',header.ti(0,1))

  c_s_to_a, istime, sstime
  c_s_to_a, ietime, setime
  c_a_to_r, istime, basetime
  c_a_to_r, ietime, endtime

; read scalor data

  if nvar_scal gt 0 then begin

    read_flat_scalor, istime, ietime, col(0:nvar_scal-1,0), stime, 	$
		      sdata, srows, filenum = id(0:nvar_scal-1,0)

  endif

; read local time data

  read_flat_scalor, istime, ietime, local_time_col, lttime, ltdata,	$
                      ltrows, filenum = local_time_id

  ltdata_b = fltarr(ltrows(0))
  ltdata_b(*) = ltdata(0,*)
  ltdata = ltdata_b

; read latitude data

  read_flat_scalor, istime, ietime, lat_col, lattime, latdata,	$
                      latrows, filenum = lat_id

  latdata_b = fltarr(latrows(0))
  latdata_b(*) = latdata(0,*)

  londata_b = fltarr(latrows(0))
  londata_b(*) = latdata(1,*)

  latdata = latdata_b
  londata = londata_b

  print, ''
  print, '1. Altitude vs. Geographic Latitude'
  print, '2. Altitude vs. Geographic Longitude'
  print, '3. Altitude vs. Magnetic Latitude'
  print, '4. Altitude vs. Geographic Londitude'
  print, '5. Geographic Latitude vs. Geographic Longitude (1 height)'
  print, '6. Geographic Latitude vs. Geographic Local Time (1 height)'
  print, '7. Magnetic Latitude vs. Magnetic Longitude (1 height)'
  print, '8. Magnetic Latitude vs. Magnetic Local Time (1 height)'
  plottype = fix(ask('type of plot','1'))

  if plottype eq 1 then begin
    xrange = [-90,90]
    xtitle_save = 'Geographic Latitude'
  endif

  if plottype eq 2 then begin
    latdata = (londata+360.0) mod 360.0
    xrange = [0,360]
    xtitle_save = 'Geographic Longitude'
  endif

  if plottype ge 3 then begin

    nlat = latrows(0)
    nhgt = n_elements(heights)

    mlat = fltarr(nlat)
    mlon = fltarr(nlat)
    mmlt = fltarr(nlat)

    if plottype ge 5 then begin

      mh = 0.0
      read, 'Enter height in km to display : ',mh
      dis        = abs(heights-mh)
      loc_height = where(dis eq min(dis))
      loc_height = loc_height(0)
      mh         = heights(loc_height)
      print, 'Closest height is ',mh,' km. Adjusting to this height.'
      h_string = 'Height : '+tostr(fix(mh))+' km'
      minlat = 0.0
      read, 'Enter minimum latitude to plot : ',minlat
      maxrange = 90.0 - minlat

      nlats = 0
      nlongs = 0

      read,'Enter number of lats, azimuths : ', nlats, nlongs

    endif else mh = mean(heights)

    if (plottype eq 3) or (plottype eq 4) or (plottype ge 7) then begin

      print, 'converting to magnetic coordinates...'

      i = long(0)

      for i=long(0),long(nlat-1) do begin
        ctime = lattime(0,i)
        c_r_to_a, itime, ctime
        year   = itime(0)+1900
        if year lt 1965 then year = year + 100
        month  = itime(1)
        day    = itime(2)
        hour   = itime(3)
        minute = itime(4)
        second = itime(5)

        utbt   = float(jday(year,month,day)-1)*24.0*3600.0 + hour*3600.0 + $
             minute*60.0 + second

        mag_pos   = cnvcoord(latdata(i),londata(i),mh)
        mag_mlt   = mlt(year, long(utbt), mag_pos(1))
        mlat(i)   = mag_pos(0)
        mlon(i)   = mag_pos(1)
        mmlt(i)   = mag_mlt
      endfor

      mlon = (mlon + 360.0) mod 360.0

    endif 

    if plottype eq 3 then begin
      xrange = [-90,90]
      xtitle_save = 'Magnetic Latitude'
      latdata(*) = mlat(*)
    endif

    if plottype eq 4 then begin
      xrange = [0,360]
      xtitle_save = 'Magnetic Longitude'
      latdata(*) = mlon(*)
    endif

    if plottype ge 5 then begin

      xrange = [0,360]

      if (plottype eq 5) then begin
        h_string = h_string + ' ; Geographic Coordinates'
        xtitle_save = 'Geographic Longitude'
      endif

      if (plottype eq 6) then begin
        h_string = h_string + ' ; Geographic Coordinates'
        xtitle_save = 'Geographic Longitude'
	londata(*) = ltdata(*)*360.0/24.0
      endif

      if (plottype eq 7) then begin
        h_string = h_string + ' ; Magnetic Coordinates'
        xtitle_save = 'Magnetic Longitude'
        latdata(*) = mlat(*)
        londata(*) = mlon(*)
      endif

      if (plottype eq 8) then begin
        h_string = h_string + ' ; Magnetic Coordinates'
        xtitle_save = 'Magnetic Longitude'
        latdata(*) = mlat(*)
        londata(*) = mmlt(*)*360.0/24.0
      endif

    endif

  endif

; read vector data

  read_flat_vector, istime, ietime, col(0:nvar_vect-1,1), vtime, vdata,	$
                      vrows, filenum = id(0:nvar_vect-1,1)

  ny = n_elements(vdata(0,0,*))

; break apart haloe data into sunrise and sunset data
;   right now we have a super simple way to determine sunrise and sunset
;   but it isn't right, so eventually we will have to come up with something
;   different

  sr_loc = where(ltdata lt 12.0, count_sr)
  ss_loc = where(ltdata gt 12.0, count_ss)

  if nvar_scal gt 0 then begin

    sdata_b = fltarr(2,nvar_scal,n_elements(sdata(0,*)))
    stime_b = dblarr(2,nvar_scal,n_elements(sdata(0,*)))

  endif

  logged = intarr(nvar_vect)

  if plottype lt 5 then begin

    nbins = (xrange(1) - xrange(0))/2.0
    if nbins gt 100 then nbins = nbins/4
    end_j = 2
    binned = fltarr(end_j+1,nvar_vect,nbins,ny)

  endif else begin

    end_j = 1
    outbin = fltarr(nlongs,nlats)
    binned = fltarr(end_j+1,nvar_vect,nlongs,nlats)
    logged(*) = 1

  endelse

  for j=0,end_j do begin

    for i=0,nvar_vect-1 do begin

      if (j eq 0) and (plottype lt 5) then 				$
        loc = where(ltdata lt 12.0, count)

      if (j eq 0) and (plottype ge 5) then 				$
        loc = where(latdata ge minlat, count)

      if (j eq 1) and (plottype lt 5) then 				$
        loc = where(ltdata gt 12.0, count)

      if (j eq 1) and (plottype ge 5) then 				$
        loc = where(-1.0*latdata ge minlat, count)

      if (j eq 2) and (plottype lt 5) then 				$
        loc = where(ltdata ge 0.0, count)

      if count gt 0 then begin

	vdata_b   = fltarr(count,ny)
	latdata_b = fltarr(count)
	londata_b = fltarr(count)

	vdata_b(0:count-1,0:ny-1) = vdata(i,loc,0:ny-1)
	latdata_b(0:count-1)      = latdata(loc)
	londata_b(0:count-1)      = londata(loc)

        if (j eq 0) and (plottype lt 5) then begin

          locnz = where(vdata_b ne 0,countnz)

          if countnz gt 0 then begin

            mm_temp = log(mm(vdata_b(locnz)))

            if ((mm_temp(1) - mm_temp(0)) lt 2.0) then begin
              print, 'Variable ',strcompress(header.na(id(i,1),col(i,1))), $
		     ' does not appear to need to be "logged".  '
              print, 'Would you like this variable to be "logged"?'
              que = ask('y or n','n')
              if strmid(que,0,1) eq 'n' then logged(i) = 1
            endif

          endif

        endif

        if plottype lt 5 then begin

          if logged(i) eq 0 then 					$
            bin_y, vdata_b, latdata(loc), nbins, outbin, 		$
	      smooth=10, range = xrange, /logbin			$
	  else								$
            bin_y, vdata_b, latdata(loc), nbins, outbin, 		$
	      smooth=10, range = xrange

        endif else begin

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
	      loc = where( (fac*latdata_b ge lat1) and			$
			   (fac*latdata_b lt lat2) and			$
			   (londata_b ge lon1) and			$
			   (londata_b lt lon2) and			$
			   (vdata_b(*,loc_height) ne 0.0),count)
	      if count gt 0 then outbin(ix,iy) = mean(vdata_b(loc,loc_height))
	    endfor
          endfor

	  y = outbin*0.0 + 1.0
	  loc = where(outbin eq 0.0,count)
	  smooth, outbin, y, 0.0
	  if count gt 0 then outbin(loc) = 0.0

        endelse

        binned(j,i,*,*) = outbin

      endif

    endfor

  endfor

; save a version of binned to calculate differences

  binned_save    = binned

; read in color table

  ct_dir = getenv('IDL_EXTRAS')
  ctname = ct_dir+'joule4.ct'
  readct,ncolors,ctname

; need to take log of the binned vector data, and scale the data :

  minmax   = fltarr(nvar_vect,2)
  minmax_d = fltarr(nvar_vect,2)

  if plottype lt 5 then begin
    temp = fltarr(end_j+1,nbins,ny)
  endif else begin
  endelse

  temp = fltarr(end_j+1,n_elements(binned(0,0,*,0)),			$
		        n_elements(binned(0,0,0,*)))

  for i = 0,nvar_vect-1 do begin

    temp(*,*,*) = binned(*,i,*,*)

    loc = where(temp ne 0,count)

    if count gt 0 then begin

      if logged(i) eq 0 then minmax(i,*) = log(mm(temp(loc)))		$
      else minmax(i,*) = mm(temp(loc))

      print, 'Variable ',header.na(id(i,1),col(i,1)),' : '
      print, '   Minimum : ',minmax(i,0)
      print, '   Maximum : ',minmax(i,1)
      read,  '   Enter minimum, maximum to scale to : ',mini,maxi
      minmax(i,*) = [mini,maxi]

    endif else minmax(i,*) = [-1,1]

    if logged(i) eq 0 then begin

      loc = where((temp lt 10.0^minmax(i,0)) and (temp ne 0.0),count)
      if count gt 0 then temp(loc) = 10.0^minmax(i,0)

      loc = where((temp gt 10.0^minmax(i,1)) or (temp eq 0.0),count)
      if count gt 0 then temp(loc) = 10.0^minmax(i,1)

      temp2 = temp

      temp = float(ncolors)*(log(temp)-minmax(i,0))/(minmax(i,1)-minmax(i,0))

    endif else begin

      loc = where((temp lt minmax(i,0)) and (temp ne 0.0),count)
      if count gt 0 then temp(loc) = minmax(i,0)

      loc = where((temp gt minmax(i,1)) or (temp eq 0.0),count)
      if count gt 0 then temp(loc) = minmax(i,1)

      temp2 = temp

      temp = float(ncolors)*(temp-minmax(i,0))/(minmax(i,1)-minmax(i,0))

    endelse

    binned(*,i,*,*) = temp

    binned_save(*,i,*,*) = temp2

  endfor

  dx = 1.25*float(!d.y_ch_size)/float(!d.y_size)

  ppp = 2
  space = dx*2.0

  pos_space,ppp,space,sizes, ny = ppp

  if (plottype lt 5) then ytitle = 'Pressure' else ytitle = 'Latitude'

  c_r_to_a, istime, basetime
  c_r_to_a, ietime, endtime
  istime([3,4,5]) = 0
  time_axis, istime, ietime, srtime, ertime,        $
             xtickname, dates, xtickvalue, xminor, xtickn

  c_a_to_s, istime, sstime
  c_a_to_s, ietime, setime
  dates = 'Dates included : '+strmid(sstime,0,9)+' to '+strmid(setime,0,9)

; plot vectors

  for var = 0,nvar_vect-1 do begin

    for j=0,end_j do begin

      if ((j mod ppp) eq 0) then begin
        plotdumb
        xyouts, 0.0, -dx*2.0, dates, /norm
        if (plottype ge 5) then xyouts, 0.0, -dx*4.0, h_string, /norm
      endif

      if (j mod ppp lt ppp-1) and (j ne 2) then begin

	xtitle = ' '
	xtickname = strarr(10)+' '

      endif else begin

	xtitle = xtitle_save
	xtickname = strarr(10)+''

      endelse

      if (j eq 0) and (plottype lt 5) then title = 'Sunrise'
      if (j eq 0) and (plottype ge 5) then title = 'Northern Hemisphere'
      if (j eq 1) and (plottype lt 5) then title = 'Sunset'
      if (j eq 1) and (plottype ge 5) then title = 'Southern Hemisphere'
      if j eq 2 then title = 'Sunrise + Sunset'

      if (id(var,1) eq 0) and (plottype lt 5)  then 			$
	yrange = [1.0e+3,1.0e-6]

      if (plottype lt 5) then						$ 
        get_position, ppp, space, sizes, j mod ppp, pos, /rect		$
      else get_position, ppp, space, sizes, j mod ppp, pos

      pos([0,2]) = pos([0,2]) - space/2.0

      if (plottype lt 5) then begin

        tv, binned(j,var,*,*), pos(0), pos(1), 			$
            xsize=pos(2)-pos(0), ysize=pos(3)-pos(1), /norm

        plot_io, xrange, yrange,		 			$
                 xtickname = xtickname,	 				$
	         xtitle = xtitle,					$
                 xstyle = 1, pos = pos, /noerase, /nodata,		$
	         ystyle = 1, yrange = yrange,				$
	         ytitle = ytitle

      endif else begin

	plot, [-maxrange,maxrange],[-maxrange,maxrange],		$
	      xstyle = 5, ystyle = 5, pos = pos,			$
	      /nodata, /noerase

	dlong = 360.0/(nlongs)
	dlat  = (90.0-minlat)/(nlats)
	if j eq 0 then fac = 1.0 else fac = -1.0

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
	    polyfill, x, y, color = binned(j,var,ix,iy)
	  endfor
        endfor

        if (plottype eq 5) or (plottype eq 7) then begin
          plotmlt, maxrange, /longs
        endif else begin
          plotmlt, maxrange
        endelse

;        plotmlt, maxrange, /no00, /no06, /no12, /no18, /white

      endelse

      xyouts, pos(0), pos(3)+0.01, title, /norm

      if con eq 2 then begin

	ma = minmax(var,1)
	mi = minmax(var,0)

	if ((ma-mi) gt 20.0) then 					$
	  dl = float(fix(fix((ma-mi)/20.0)/5.0))*5.0+5.0		$
	else								$
	  dl = float(fix(fix((ma-mi)/20.0)/1.0))*1.0+1.0

	mi = float(fix(mi/dl))*dl - dl

	levels = findgen(30)*dl + mi

        contour, binned_save(j,var,*,*), levels = levels, 		$
		 min_value = minmax(var,0)*1.01,			$
		 max_value = minmax(var,1)*0.99,			$
		 pos = pos, /noerase, /follow,				$
		 xstyle = 5, ystyle = 5

      endif

; plot color table

      ct_title = header.na(id(var,1),col(var,1))
      posct = [0.90,pos(1),0.92,pos(3)]
      plotct, ncolors, posct, minmax(var,*), ct_title, /right

    endfor

  endfor

;  if (!d.name eq 'PS') then begin

;    device, /close
;    psfile2 = strmid(psfile,0,strpos(psfile,'.ps')) + '_bw.ps'
;    setdevice, psfile2,'l',4,0.95

;  endif

endif

if (!d.name eq 'PS') then begin

  device, /close
  set_plot, 'X'

endif

end
