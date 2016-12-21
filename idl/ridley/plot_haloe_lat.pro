
common ffinfo, header

; ----------------------------------------------------------------
; BIG ASSUMPTIONS :

local_time_id  = [0]
local_time_col = [2]

lat_id  = [0]
lat_col = [0]

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
  latdata = latdata_b

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

  xrange = [-75,75]
  xtitle = 'Latitude'

  logged = intarr(nvar_vect)

  nbins = (xrange(1) - xrange(0))/2.0

  binned = fltarr(3,nvar_vect,nbins,ny)

  for j=0,2 do begin

    for i=0,nvar_vect-1 do begin

      if j eq 0 then 							$
        loc = where(ltdata lt 12.0, count)

      if j eq 1 then 							$
        loc = where(ltdata gt 12.0, count)

      if j eq 2 then 							$
        loc = where(ltdata ge 0.0, count)

      if count gt 0 then begin

	vdata_b = fltarr(count,ny)
	vdata_b(0:count-1,0:ny-1) = vdata(i,loc,0:ny-1)

        if j eq 0 then begin

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

        if logged(i) eq 0 then 						$
          bin_y, vdata_b, latdata(loc), nbins, outbin, 			$
	     smooth=10, range = xrange, /logbin				$
	else								$
          bin_y, vdata_b, latdata(loc), nbins, outbin, 			$
	     smooth=10, range = xrange

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
  temp = fltarr(3,nbins,ny)

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

  ppp = 3
  space = dx*2.0

  pos_space,ppp,space,sizes, ny = ppp

  ytitle = 'Pressure'

  c_r_to_a, istime, basetime
  c_r_to_a, ietime, endtime
  istime([3,4,5]) = 0
  time_axis, istime, ietime, srtime, ertime,        $
             xtickname, dates, xtickvalue, xminor, xtickn

  dates = 'Dates included : '+strmid(dates,0,strpos(dates,'U')-1)

; plot vectors

  for var = 0,nvar_vect-1 do begin

    for j=0,2 do begin

      if j eq 0 then plotdumb

      if j lt 2 then begin

	xtitle = ' '
	xtickname = strarr(10)+' '

      endif else begin

	xtitle = 'Latitude'
	xtickname = strarr(10)+''

      endelse

      if j eq 0 then title = 'Sunrise'
      if j eq 1 then title = 'Sunset'
      if j eq 2 then title = 'Sunrise + Sunset'

      if id(var,1) eq 0 then yrange = [1.0e+3,1.0e-6]			$
      else yrange = [0.0,1.0]

      get_position, ppp,space, sizes, j, pos, /rect
      pos([0,2]) = pos([0,2]) - space/2.0

      tv, binned(j,var,*,*), pos(0), pos(1), 			$
          xsize=pos(2)-pos(0), ysize=pos(3)-pos(1), /norm

      plot_io, xrange, yrange,		 			$
            xtickname = xtickname,	 				$
	    xtitle = xtitle,						$
            xstyle = 1, pos = pos, /noerase, /nodata,			$
	    ystyle = 1, yrange = yrange,				$
	    ytitle = ytitle

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
      posct = [0.98,pos(1),1.0,pos(3)]
      plotct, ncolors, posct, minmax(var,*), ct_title, /right

    endfor

    xyouts, 0.0, -dx*2.0, dates, /norm

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
