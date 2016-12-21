
common ffinfo, header

; ----------------------------------------------------------------
; BIG ASSUMPTIONS :

local_time_id  = [0]
local_time_col = [2]

lon_id  = [0]
lon_col = [1]

lat_id  = [0]
lat_col = [0]

; ----------------------------------------------------------------

; Pressure Levels :

pressure = [1.000e+03, 6.813e+02, 4.642e+02, 3.162e+02, 2.154e+02, 1.468e+02, $
            1.000e+02, 6.813e+01, 4.642e+01, 3.162e+01, 2.154e+01, 1.468e+01, $
            1.000e+01, 6.813e+00, 4.642e+00, 3.162e+00, 2.154e+00, 1.468e+00, $
            1.000e+00, 6.813e-01, 4.642e-01, 3.162e-01, 2.154e-01, 1.468e-01, $
            1.000e-01, 6.813e-02, 4.642e-02, 3.162e-02, 2.154e-02, 1.468e-02, $
            1.000e-02, 6.813e-03, 4.642e-03, 3.162e-03, 2.154e-03, 1.468e-03, $
            1.000e-03, 6.813e-04, 4.642e-04, 3.162e-04, 2.154e-04, 1.468e-04, $
            1.000e-04, 6.813e-05, 4.642e-05, 3.162e-05, 2.154e-05, 1.468e-05, $
            1.000e-05, 6.813e-06, 4.642e-06, 3.162e-06, 2.154e-06, 1.468e-06, $
            1.000e-06]

; Height Levels :

heights =  [1.102e-01, 9.127e-01, 6.125e+00, 8.822e+00, 1.133e+01, 1.378e+01, $
            1.622e+01, 1.867e+01, 2.112e+01, 2.360e+01, 2.610e+01, 2.864e+01, $
            3.121e+01, 3.382e+01, 3.651e+01, 3.929e+01, 4.216e+01, 4.512e+01, $
            4.818e+01, 5.127e+01, 5.432e+01, 5.727e+01, 6.014e+01, 6.292e+01, $
            6.562e+01, 6.823e+01, 7.077e+01, 7.323e+01, 7.564e+01, 7.800e+01, $
            8.030e+01, 8.256e+01, 8.477e+01, 8.693e+01, 8.910e+01, 9.126e+01, $
            9.344e+01, 9.562e+01, 9.783e+01, 1.001e+02, 1.024e+02, 1.049e+02, $
            1.075e+02, 1.103e+02, 1.136e+02, 1.175e+02, 1.221e+02, 1.275e+02, $
            1.338e+02, 1.411e+02, 1.495e+02, 1.590e+02, 1.695e+02, 1.813e+02, $
            1.941e+02]

haloefile = ''
read, 'Enter Flatfile name for HALOE data : ', haloefile

paramfile = ''
read, 'Enter Flatfile name for index data (return for none) : ', paramfile

psfile = ''
read, 'Enter ps filename (return for screen) : ',psfile

if strlen(psfile) gt 0 then setdevice, psfile,'l',4,0.9			$
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

  ndays = fix((endtime-basetime)/(24.0*3600.0) + 0.5)

  mh = 0.0
  read, 'Enter height in km to display : ',mh
  dis        = abs(heights-mh)
  loc_height = where(dis eq min(dis))
  lh         = loc_height(0)
  mh         = heights(lh)
  print, 'Closest height is ',mh,' km. Adjusting to this height.'
  h_string = 'Height : '+tostr(fix(mh))+' km'

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
  latdata = latdata_b

; read longitude data

  read_flat_scalor, istime, ietime, lon_col, lontime, londata,	$
                      lonrows, filenum = lon_id

  londata_b = fltarr(lonrows(0))
  londata_b(*) = londata(0,*)
  londata = londata_b

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

  if lon_col(0) eq 1 then begin
    xrange = [0,360]
    xtitle = 'Longitude'
    yrange = [30,60]
    ytitle_save = 'Latitude'
  endif

  nbins   = (xrange(1) - xrange(0))/30
  nbins_y = (yrange(1) - yrange(0))/5

  dx = (xrange(1)-xrange(0))/nbins
  dy = (yrange(1)-yrange(0))/nbins_y

  binned = fltarr(2,nvar_vect,ndays,nbins,nbins_y)

  for j=0,1 do begin

    for k=0,ndays do for i=0,nvar_vect-1 do begin

      st = basetime + double(k)*24.0*3600.0
      et = st + 24.0*3600.0

      if j eq 0 then 							$
        loc = where(ltdata lt 12.0 and 					$
                    latdata ge yrange(0) and				$
                    latdata lt yrange(1) and				$
                    vtime(i,*) ge st and 				$
		    vtime(i,*) lt et, count)

      if j eq 1 then 							$
        loc = where(ltdata gt 12.0 and 					$
                    latdata ge yrange(0) and				$
                    latdata lt yrange(1) and				$
                    vtime(i,*) ge st and 				$
		    vtime(i,*) lt et, count)

      if count gt 0 then begin

	vdata_b = fltarr(count,ny)
	vdata_b(0:count-1,0:ny-1) = vdata(i,loc,0:ny-1)

        londata(loc) = (londata(loc) mod 360.0)

        outbin       = fltarr(nbins,nbins_y)
        outbin_count = intarr(nbins,nbins_y)

	for pts = 0,count-1 do begin
          outbin(londata(loc(pts))/dx,					$
                 (latdata(loc(pts))-yrange(0))/dy) = vdata_b(pts,lh)
          outbin_count(londata(loc(pts))/dx,				$
                      (latdata(loc(pts))-yrange(0))/dy) = 		$
                      outbin_count(londata(loc(pts))/dx,		$
                      (latdata(loc(pts))-yrange(0))/dy) + 1
	endfor

	loc = where(outbin_count gt 0, count)

	if count gt 0 then outbin(loc) = outbin(loc)/float(outbin_count(loc))

;        bin_y, vdata_b, londata(loc), nbins, outbin, 			$
;	     smooth=10, range = xrange, /logbin,			$
;	     /rotate

        binned(j,i,k,*,*) = outbin

      endif

    endfor

  endfor

; save a version of binned to calculate differences

  binned_save    = binned
  binned_save_2  = binned

; read in color table

  ct_dir = getenv('IDL_EXTRAS')
  ctname = ct_dir+'joule4.ct'
  readct,ncolors,ctname

; need to take log of the binned vector data, and scale the data :

  minmax   = fltarr(nvar_vect,2)
  minmax_d = fltarr(nvar_vect,2)
  temp = fltarr(2,ndays,nbins,nbins_y)

  for i = 0,nvar_vect-1 do begin

    temp(*,*,*,*) = binned(*,i,*,*,*)

    loc = where(temp ne 0,count)

    if count gt 0 then begin

      minmax(i,*) = log(mm(temp(loc)))

      print, 'Variable ',header.na(id(i,1),col(i,1)),' : '
      print, '   Minimum : ',minmax(i,0)
      print, '   Maximum : ',minmax(i,1)
      read,  '   Enter minimum, maximum to scale to : ',mini,maxi
      minmax(i,*) = [mini,maxi]

    endif else minmax(i,*) = [-1,1]

    loc = where((temp lt 10.0^minmax(i,0)) and (temp ne 0.0),count)
    if count gt 0 then temp(loc) = 10.0^minmax(i,0)

    loc = where((temp gt 10.0^minmax(i,1)) or (temp eq 0.0),count)
    if count gt 0 then temp(loc) = 10.0^minmax(i,1)

    temp = float(ncolors)*(log(temp)-minmax(i,0))/(minmax(i,1)-minmax(i,0))

    binned(*,i,*,*,*) = temp

  endfor

; read in pem ionization rates

  list = findfile("/d/ridley/d.data/d.haloe/d.flat/d.energy/energy*.ascii")
  nfile = n_elements(list)

print, nfile

  pembin   = intarr(nfile,nbins,nbins_y)
  pemcount = intarr(nfile,nbins,nbins_y)
  pemalt      = strarr(nfile)
  pemday      = intarr(nfile)

  line = ''

  for i=0,nfile-1 do begin

    openr,1,list(i)
    readf,1,line
    pemday(i) = fix(strmid(line,5,3))
    readf,1,line
    pemalt(i) = strmid(line,11,8)

    while not eof(1) do begin
      readf,1,line
      lon = float(strmid(line,6,6))+2.5
      lat = float(strmid(line,32,6))+2.5
      if lat ge yrange(0) and lat lt yrange(1) then begin
        pembin(i,lon/dx,(lat-yrange(0))/dy) = float(strmid(line,60,11))
        pemcount(i,lon/dx,(lat-yrange(0))/dy) = 			$
          pemcount(i,lon/dx,(lat-yrange(0))/dy) + 1
      endif
    endwhile

    close,1

  endfor  

  loc = where(pemcount gt 0 and pembin gt 0)
  locz = where(pemcount eq 0 or pembin le 0)
  pembin(loc) = pembin(loc)/pemcount(loc)
  pembin(loc) = ncolors*(alog10(pembin(loc))+3)/8.0
  pembin(locz) = ncolors

  ppp = 4
  space = 0.075

  pos_space,ppp,space,sizes, ny = ppp/2.0

  for j=0,0 do begin

; plot vectors

    for var = 0,nvar_vect-1 do begin

      for k=0,ndays-1,2 do begin

	if k eq 0 then ik = 0
	if k eq 2 then ik = 1

        st = basetime + double(k)*24.0*3600.0
	c_r_to_a, itime, st
	c_a_to_s, itime, title
	title = strmid(title,0,9)

        plot_num = ik*2 mod ppp

	ytitle = ytitle_save
	ytickname = strarr(10)+''

        if (plot_num eq ppp-1) or (k eq ndays-1) then begin
          xtickname2 = strarr(10)+''
          xtitle2 = xtitle
        endif else begin
          xtickname2 = strarr(10)+' '
          xtitle2 = ' '
        endelse

        get_position, ppp,space, sizes, plot_num, pos, /rect
	pos([0,2]) = pos([0,2]) - space/2.0

        if plot_num eq 0 then begin

          plotdumb

          if j eq 0 then xyouts, pos(2), pos(3)+0.005, 			$
            'SunRise ('+h_string+')', /norm, alignment = 1.0
          if j eq 1 then xyouts, pos(2), pos(3)+0.005, 			$
            'SunSet ('+h_string+')', /norm, alignment = 1.0

        endif

	if (k lt ndays) then begin

          tv, binned(j,var,k,*,*), pos(0), pos(1), 			$
              xsize=pos(2)-pos(0), ysize=pos(3)-pos(1), /norm

          plot, xrange, yrange,			 			$
              xtickname = xtickname2,	 				$
              xstyle = 1, pos = pos, /noerase, /nodata,			$
	      ystyle = 1, yrange = yrange,				$
	      ytitle = ytitle,						$
	      ytickname = ytickname,					$
              xtitle = xtitle2

          xyouts, pos(0), pos(3)+0.005, title, /norm, alignment = 0.0

	endif

; plot color table 

        ct_title = header.na(id(var,1),col(var,1))
        posct = [pos(2)+0.005,pos(1),pos(2)+0.015,pos(3)]
        plotct, ncolors, posct, minmax(var,*), ct_title, /right

; plot ionization rates

	if k eq 0 then kk = 1
	if k eq 2 then kk = 2

	itime(0) = 1993
        itime(1) = 1
        itime(2) = pemday(kk)
	c_a_to_r, itime, temp
        c_r_to_a, itime, temp
	c_a_to_s, itime, title
	title = strmid(title,0,9)

        plot_num = (ik*2+1) mod ppp

	ytitle = ' '
	ytickname = strarr(10)+' '

        if (plot_num eq ppp-1) or (k eq ndays-1) then begin
          xtickname2 = strarr(10)+''
          xtitle2 = xtitle
        endif else begin
          xtickname2 = strarr(10)+' '
          xtitle2 = ' '
        endelse

        get_position, ppp,space, sizes, plot_num, pos, /rect
;	pos([0,2]) = pos([0,2]) - space/2.0

        xyouts, pos(2), pos(3)+0.005, 			$
          'Ionization ('+pemalt(kk)+')', /norm, alignment = 1.0

	if (k lt ndays) then begin

          tv, pembin(kk,*,*), pos(0), pos(1), 			$
              xsize=pos(2)-pos(0), ysize=pos(3)-pos(1), /norm

          plot, xrange, yrange,			 			$
              xtickname = xtickname2,	 				$
              xstyle = 1, pos = pos, /noerase, /nodata,			$
	      ystyle = 1, yrange = yrange,				$
	      ytitle = ytitle,						$
	      ytickname = ytickname,					$
              xtitle = xtitle2

          xyouts, pos(0), pos(3)+0.005, title, /norm, alignment = 0.0

	endif

; plot color table 

        ct_title = 'Log Ionizations/(cm3 s)'
        posct = [pos(2)+0.005,pos(1),pos(2)+0.015,pos(3)]
        plotct, ncolors, posct, [-3.0,5.0], ct_title, /right

      endfor

    endfor

  endfor

  if (!d.name eq 'PS') then begin

    device, /close
    psfile2 = strmid(psfile,0,strpos(psfile,'.ps')) + '_bw.ps'
    setdevice, psfile2,'l',4,0.95

  endif

  c_r_to_a, istime, basetime
  c_r_to_a, ietime, endtime
  istime([3,4,5]) = 0
  time_axis, istime, ietime, srtime, ertime,        $
             xtickname, xtitle_bw, xtickvalue, xminor, xtickn
  times = findgen(ndays)*24.0*3600.0 + 12.0*3600.0
  dx = 1.25*float(!d.y_ch_size)/float(!d.y_size)

;  bins = [6,24]
;  alts = [45,45]

  for i = 0,n_elements(bins)-1 do begin

    bin = bins(i)
    alt = alts(i)

    pressure = float(alt)*(log(yrange(1)) - log(yrange(0)))/float(ny) + $
               log(yrange(0))
    pressure = 'Pressure : '+tostrf(10.0^pressure)

    ll = float(bin)*(xrange(1) - xrange(0))/float(nbins) + xrange(0)
    ll = xtitle + ' : ' + tostr(ll)

    for j=0,1 do begin

      if j eq 0 then sun = 'Sunrise Data'
      if j eq 1 then sun = 'Sunset Data'

      for var = 0,nvar_vect-1 do begin

        variable = 'Variable : '+header.na(id(var,1),col(var,1))

        pos = [0.0,0.01,0.99,1.0-4.5*dx]

        plot, times, binned_save_2(j,var,*,bin,alt), 			$
              xtickname = xtickname,	 				$
              xtickv = xtickvalue, xticks = xtickn, 			$
              xminor = xminor, xtitle=xtitle_bw,			$
	      xrange = [srtime, ertime],			        $
	      ytitle = 'Percent Difference', pos = pos

        xyouts, 0.02, 1.0-1.0*dx, sun, /norm
        xyouts, 0.02, 1.0-2.0*dx, variable, /norm
        xyouts, 0.02, 1.0-3.0*dx, pressure, /norm
        xyouts, 0.02, 1.0-4.0*dx, ll, /norm

      endfor

    endfor

  endfor

endif

if (!d.name eq 'PS') then begin

  device, /close
  set_plot, 'X'

endif

end
