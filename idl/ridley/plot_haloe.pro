
common ffinfo, header

; ----------------------------------------------------------------
; BIG ASSUMPTIONS :

local_time_id  = [0]
local_time_col = [2]

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

datafile = strarr(6)
datafile(0) = haloefile

done = 0
nfile = 1

while not done do begin

  que = ''
  read, 'Enter Flatfile name for index data (return to end) : ', que

  if strlen(que) gt 0 then begin
    datafile(nfile) = que
    nfile = nfile + 1
  endif else done = 1

endwhile

psfile = ''
read, 'Enter ps filename (return for screen) : ',psfile

if strlen(psfile) gt 0 then setdevice, psfile,'p',4,0.85		$
else window,0

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

  if nvar_scal gt 0 then begin

; read scalor data

    read_flat_scalor, istime, ietime, col(0:nvar_scal-1,0), stime, 	$
		      sdata, srows, filenum = id(0:nvar_scal-1,0)

  endif

; read local time data

  read_flat_scalor, istime, ietime, local_time_col, lttime, ltdata,	$
                      ltrows, filenum = local_time_id

  if nvar_vect gt 0 then begin

; read vector data

    read_flat_vector, istime, ietime, col(0:nvar_vect-1,1), vtime, 	$
		      vdata, vrows, filenum = id(0:nvar_vect-1,1)

  endif

  ny = n_elements(vdata(0,0,*))

; break apart haloe data into sunrise and sunset data
;   right now we have a super simple way to determine sunrise and sunset
;   but it isn't right, so eventually we will have to come up with something
;   different

  sr_loc = where(ltdata lt 12.0, count_sr)
  ss_loc = where(ltdata gt 12.0, count_ss)

  nbins = max([count_sr,count_ss])

  if nvar_scal gt 0 then begin
    sdata_b = fltarr(2,nvar_scal,n_elements(sdata(0,*)))
    stime_b = dblarr(2,nvar_scal,n_elements(sdata(0,*)))
  endif

  if nvar_vect gt 0 then begin
    vdata_b = fltarr(n_elements(vdata(0,*,0)), ny)
    vtime_b = dblarr(n_elements(vdata(0,*,0)))
    binned = fltarr(2,nvar_vect,nbins,ny)
    logged = intarr(nvar_vect)
  endif

  range = [basetime, endtime]

  for j=0,1 do begin

    if j eq 0 then begin
      loc = sr_loc
      count = count_sr
    endif

    if j eq 1 then begin
      loc = ss_loc
      count = count_ss
    endif

    if (count gt 0) then begin

      for i=0,nvar_scal-1 do begin

        if id(i,0) eq 0 then begin
          sdata_b(j,i,0:count-1) = sdata(i,loc)
          stime_b(j,i,0:count-1) = stime(i,loc)
        endif else begin
          sdata_b(j,i,*) = sdata(i,*)
          stime_b(j,i,*) = stime(i,*)
        endelse

      endfor

      for i=0,nvar_vect-1 do begin

        if id(i,1) eq 0 then begin
          vdata_b(0:count-1,0:ny-1) = vdata(i,loc,0:ny-1)
          vtime_b(0:count-1)        = vtime(i,loc)
        endif else begin
          vdata_b(*,0:ny-1)         = vdata(i,*,0:ny-1)
          vtime_b(*)                = vtime(i,*)
        endelse

        if j eq 0 then begin

          locnz = where(vdata_b ne 0,countnz)

          if countnz gt 0 then begin

            mm_temp = log(mm(vdata_b(locnz)))

            if ((mm_temp(1) - mm_temp(0)) lt 4.0) then begin
              print, 'Variable ',strcompress(header.na(id(i,1),col(i,1))), $
		     ' does not appear to need to be "logged".  '
              print, 'Would you like this variable to be "logged"?'
              que = ask('y or n','n')
              if strmid(que,0,1) eq 'n' then logged(i) = 1
            endif

          endif

        endif

        if logged(i) eq 0 then 						$
          bin_y, vdata_b, vtime_b, nbins, outbin, 			$
	         smooth=10, range = range, /logbin			$
        else								$
          bin_y, vdata_b, vtime_b, nbins, outbin, 			$
	         smooth=10, range = range

        binned(j,i,*,*) = outbin

      endfor

    endif

  endfor

; determine what type of plot we want:

  if nvar_vect gt 0 then begin
    take_sd = fltarr(nvar_vect)
    take_mn = fltarr(nvar_vect)
  endif

  for i = 0,nvar_vect-1 do begin

    print, 'Variable ',header.na(id(i,1),col(i,1)),' : '
    print, '  1. Normal'
    print, '  2. Averaged (log)'
    print, '  3. Averaged (no log)'
    print, '  4. Standard Deviation'
    print, '  5. Standard Deviation (as percent of mean)'
    plottype = fix(ask('type of plot','1'))

    if plottype eq 2 then take_mn(i) = 1
    if plottype eq 3 then take_mn(i) = 2
    if plottype eq 4 then take_sd(i) = 1
    if plottype eq 5 then take_sd(i) = 2

  endfor

  ndays = fix((endtime-basetime)/(24.0*3600.0) + 0.5)

; read in color table

  ct_dir = getenv('IDL_EXTRAS')
  ctname = ct_dir+'joule4.ct'
  readct,ncolors,ctname

; need to take log of the binned vector data, and scale the data :

  temp = fltarr(nbins,ny)

  if nvar_vect gt 0 then begin
    minmax = fltarr(nvar_vect,2)
    binned_save    = binned
    binned_save_mn = binned
  endif

  for j=0,1 do for i = 0,nvar_vect-1 do begin

    temp(*,*) = binned(j,i,*,*)

    nx = n_elements(temp(*,0))
    ny = n_elements(temp(0,*))

    if take_sd(i) gt 0 then begin

      std = temp*0.0
      mn = temp*0.0
      mnl = temp*0.0

      dx = nbins/ndays

      for x=0,nx-1 do begin

	sx = x-dx/2
	if sx lt 0 then sx = 0
	ex = sx + dx
	if ex gt nx-1 then begin
          ex = nx-1
	  sx = max([0,ex-dx])
	endif

	for y=0,ny-1 do begin

	  std(x,y) = stdev(temp(sx:ex,y),0.0)
          if take_sd(i) eq 2 then mn(x,y) = mean(temp(sx:ex,y),0.0)

	endfor

      endfor

      if take_sd(i) eq 2 then begin
        loc = where(mn ne 0.0,count)
        if count gt 0 then std(loc) = 100.0*std(loc)/mn(loc)
        logged(i) = 1
      endif

      temp = std

    endif

    if take_mn(i) gt 0 then begin

      mn = temp*0.0

      dx = nbins/ndays

      for x=0,nx-1 do begin

	sx = x-dx/2
	if sx lt 0 then sx = 0
	ex = sx + dx
	if ex gt nx-1 then begin
          ex = nx-1
	  sx = max([0,ex-dx])
	endif

	if take_mn(i) eq 1 then 					$
	  for y=0,ny-1 do mn(x,y) = mean(temp(sx:ex,y),0.0,/log)

	if take_mn(i) eq 2 then 					$
	  for y=0,ny-1 do mn(x,y) = mean(temp(sx:ex,y),0.0)

      endfor

      temp = mn

    endif

    loc = where(temp ne 0,count)

    if j eq 0 then begin

      if logged(i) eq 0 then minmax(i,*) = log(mm(temp(loc)))		$
      else minmax(i,*) = mm(temp(loc))

      print, 'Variable ',header.na(id(i,1),col(i,1)),' : '
      print, '   Minimum : ',minmax(i,0)
      print, '   Maximum : ',minmax(i,1)
      read,  '   Enter minimum, maximum to scale to : ',mini,maxi
      minmax(i,*) = [mini,maxi]

    endif

    binned_save_mn(j,i,*,*) = temp

    if logged(i) eq 0 then begin

      loc = where((temp lt 10.0^minmax(i,0)) and (temp ne 0.0),count)
      if count gt 0 then temp(loc) = 10.0^minmax(i,0)

      loc = where((temp gt 10.0^minmax(i,1)) or (temp eq 0.0),count)
      if count gt 0 then temp(loc) = 10.0^minmax(i,1)

      temp = float(ncolors)*(log(temp)-minmax(i,0))/(minmax(i,1)-minmax(i,0))

    endif else begin

      loc = where((temp lt minmax(i,0)) and (temp ne 0.0),count)
      if count gt 0 then temp(loc) = minmax(i,0)

      loc = where((temp gt minmax(i,1)) or (temp eq 0.0),count)
      if count gt 0 then temp(loc) = minmax(i,1)

      temp = float(ncolors)*(temp-minmax(i,0))/(minmax(i,1)-minmax(i,0))

    endelse

    binned(j,i,*,*) = temp

  endfor

  que=''
  read, 'Would you like long tick marks (y/n)? [n] ', que

  if strmid(que,0,1) eq 'y' then longticks = 1 else longticks = 0

  que=''
  read, 'Would you like sunrise and sunset data on the same page '+	$
        '(y/n)? [n] ', que

  done = 0
  pres = [-1]

  while not done do begin

    bwque = ''
    read, 'Enter pressure for black and white plots (return to end) : ',bwque
    if strlen(bwque) gt 0 then begin
      if pres(0) eq -1 then pres = float(bwque)			$
      else pres = [pres, float(bwque)]
    endif else done = 1

  endwhile

  lower = '1.0e+3'
  upper = '1.0e-6'

  rlower = float(ask('bottom pressure level',lower))
  rupper = float(ask('top pressure level',upper))

  loc_press = where(pressure le rlower and pressure ge rupper, count)

  if count eq 0 then loc_press = where(pressure ne -1.0e32)

  nvect = nvar_vect

  if strmid(que,0,1) eq 'y' then nvect = nvect*2

  nv = nvar_scal + nvect

  xmargin = 0.1
  if (nvect le 2) and (nv lt 5) then 					$
    ymargin = float(3-nvect)*0.25 					$
  else ymargin = 0.0

  ppp = nv*2
  space = 0.01

  pos_save = fltarr(nv,4)

  pos_space,ppp,space,sizes, ny = ppp

  plot_num = 0

  for i=0,nvar_scal-1 do begin

    get_position, ppp, space, sizes, i, pos, /rect, 		$
		  xmargin = xmargin, ymargin=ymargin

    pos(0) = pos(0) - xmargin + space*3.0
    pos(2) = pos(2) - xmargin + space*3.0

    pos_save(plot_num,*) = pos(*)
    plot_num = plot_num + 1

  endfor

  if nvar_scal gt 0 then 						$
    ytop = (1.0-pos(1)) + space/2.0					$
  else ytop = 0.0

  ymargin = ymargin + ytop

  ppp = nvect
  pos_space,ppp,space,sizes, ny = ppp

  for i=0,nvect-1 do begin

    get_position, ppp, space, sizes, i, pos, /rect, 		$
		  xmargin = xmargin, ymargin=ymargin

    pos(0) = pos(0) - xmargin + space*3.0
    pos(2) = pos(2) - xmargin + space*3.0
    pos(1) = pos(1) - ytop
    pos(3) = pos(3) - ytop

    pos_save(plot_num,*) = pos(*)
    plot_num = plot_num + 1

  endfor

  time_axis, istime, ietime, srtime, ertime,        $
             xtickname, xtitle, xtickvalue, xminor, xtickn
  xtickname2 = strarr(10)+' '
  xtitle2 = ' '

  plot_num = 0

  for j=0,1 do begin

    if (nv eq nvar) or (j eq 0) then begin
      plot_num = 0
      plotdumb
    endif

    if (nv eq nvar) then begin
      if j eq 0 then xyouts, 0.5, 1.0, 'SunRise', 			$
        /norm, alignment = 0.5						$
      else xyouts, 0.5, 1.0, 'SunSet', /norm, alignment = 0.5
    endif

; plot scalors

    if (nv eq nvar) or (j eq 0) then begin

      for var = 0,nvar_scal-1 do begin

        if var eq nvar-1 then begin
          xtickname_temp = xtickname
          xtitle_temp = xtitle
        endif else begin
          xtickname_temp = xtickname2
          xtitle_temp = xtitle2
        endelse

        ytitle = strcompress(header.na(id(var,0),col(var,0)),/remove_all)

        pos(*) = pos_save(plot_num,*)
        plot_num = plot_num + 1

        loc = where((stime_b(j,var,*) gt 0) and 			$
	            (sdata_b(j,var,*) ne -1.0e32), count)

        if count gt 0 then begin

          if (id(var,0) eq 0) and (nv ne nvar) then begin
            mm1 = mm(sdata_b(j,var,loc))
            loc2 = where(stime_b(1-j,var,*) gt 0, count2)
            if count2 gt 0 then mm2 = mm(sdata_b(1-j,var,loc2))		$
            else mm2 = mm1
	    yrange = mm([mm1,mm2])
	  endif else							$
            yrange = mm(sdata_b(j,var,loc))

          plot, stime_b(j,var,loc)-basetime, sdata_b(j,var,loc), 	$
              xtickname = xtickname_temp, 				$
              xtickv = xtickvalue, xticks = xtickn, 			$
              xminor = xminor, xtitle=xtitle_temp,			$
	      xrange = [srtime, ertime],		        	$
              xstyle = 1, pos = pos, /noerase,				$
	      ytitle = ytitle,						$
	      yrange = yrange

          if (id(var,0) eq 0) and (nv ne nvar) then begin

            if count2 gt 0 then 					$
	      oplot, stime_b(1-j,var,loc2)-basetime,			$
		   sdata_b(1-j,var,loc2),				$
		   linestyle = 2
	    plot, [0,1], pos = pos, /nodata, /noerase, 			$
	        xstyle=5, ystyle=5

            oplot, [0.05,0.15],[0.05,0.05]
	    xyouts, 0.16, 0.05, 'Sunrise', charsize=0.8
            oplot, [0.45,0.55],[0.05,0.05], linestyle = 2
	    xyouts, 0.56, 0.05, 'Sunset', charsize=0.8

          endif

        endif

      endfor

    endif

; plot vectors

    for var = 0,nvar_vect-1 do begin

      pos(*) = pos_save(plot_num,*)
      plot_num = plot_num + 1

      if plot_num eq nv then begin
        xtickname_temp = xtickname
        xtitle_temp = xtitle
      endif else begin
        xtickname_temp = xtickname2
        xtitle_temp = xtitle2
      endelse

      ytitle = strcompress(header.na(id(var,1),col(var,1)),/remove_all)

      if take_sd(var) then ytitle = ytitle + ' (Std. Dev.)'
      if take_mn(var) then ytitle = ytitle + ' (Mean)'

      if id(var,1) eq 0 then yrange = [rlower,rupper]			$
      else yrange = [0.0,1.0]

      tv, binned(j,var,*,loc_press), pos(0), pos(1), 			$
             xsize=pos(2)-pos(0), ysize=pos(3)-pos(1), /norm

      plot_io, [srtime,ertime], yrange,	 				$
              xtickname = xtickname_temp, 				$
              xtickv = xtickvalue, xticks = xtickn, 			$
              xminor = xminor, xtitle=xtitle_temp,			$
	      xrange = [srtime, ertime],		        	$
              xstyle = 1, pos = pos, /noerase, /nodata,			$
	      ystyle = 1, yrange = yrange,				$
	      ytitle = 'Pressure',ytickformat='(e7.1)'

      if id(var,1) eq 0 then begin
        xm = ertime + (ertime - srtime)/200.0

        if 1.3e2 lt rlower and 1.3e2 gt rupper then			$
          xyouts, xm, 1.3e2, '16', charsize = 0.8
        if 1.3e0 lt rlower and 1.3e0 gt rupper then			$
          xyouts, xm, 1.3e0, '48', charsize = 0.8
        if 1.3e-2 lt rlower and 1.3e-2 gt rupper then			$
          xyouts, xm, 1.3e-2, '80', charsize = 0.8
        if 1.3e-4 lt rlower and 1.3e-4 gt rupper then			$
          xyouts, xm, 1.3e-4, '105', charsize = 0.8
        if 1.3e-6 lt rlower and 1.3e-6 gt rupper then			$
          xyouts, xm, 1.3e-6, '220', charsize = 0.8

        if longticks then begin
	  oplot, [srtime,ertime],[1.0e2,1.0e2]
	  oplot, [srtime,ertime],[1.0e0,1.0e0]
	  oplot, [srtime,ertime],[1.0e-2,1.0e-2]
	  oplot, [srtime,ertime],[1.0e-4,1.0e-4]
        endif

        if (nv ne nvar) then begin

          plot, [0,1], pos = pos, /nodata, /noerase, 			$
	        xstyle=5, ystyle=5

	  if j eq 0 then xyouts, 0.01, 0.01, 'Sunrise', charsize=0.8	$
	  else xyouts, 0.01, 0.01, 'Sunset', charsize=0.8

        endif

      endif

      posct = [0.98,pos(1),1.0,pos(3)]

      plotct, ncolors, posct, minmax(var,*), ytitle, /right

    endfor

    if !d.name eq 'X' then prompt_for_next

  endfor

  if n_elements(pres) eq 0 then pres = [-1]

  if pres(0) ne -1 then begin

    if (!d.name eq 'PS') then begin

      device, /close
      psfile2 = strmid(psfile,0,strpos(psfile,'.ps')) + '_bw.ps'
      setdevice, psfile2,'p',4,0.95

    endif

    alts = ny*(log(pres)-log(1.0e3))/(log(1.0e-6)-log(1.0e3))

    ppp = n_elements(alts)+1
    space = 0.01

    pos_space,ppp,space,sizes, ny = ppp

    time = (findgen(nbins)+0.5)*ertime/float(nbins)

    ymargin = 0.0

    for j=0,1 do begin

; plot vectors

      for var = 0,nvar_vect-1 do begin

        if id(var,1) eq 0 then yrange = [1.0e+3,1.0e-6]			$
        else yrange = [0.0,1.0]

	plot_num = 0

        for i = 0,n_elements(alts)-1 do begin

          if plot_num eq 0 then begin

            plotdumb
            if j eq 0 then xyouts, 0.5, 1.0, 'SunRise', /norm, alignment = 0.5
            if j eq 1 then xyouts, 0.5, 1.0, 'SunSet', /norm, alignment = 0.5

	    loc = where(id(*,0) gt 0, count)

	    if count gt 0 then begin

              get_position, ppp,space, sizes, plot_num, pos, /rect, 	$
		    xmargin = xmargin, ymargin=ymargin

              pos(0) = pos(0) - xmargin + space*3.0
              pos(2) = pos(2) - xmargin + space*3.0
              xtickname_temp = xtickname2
              xtitle_temp = xtitle2

	      for ind=0,count-1 do begin

		n = loc(ind)

                ytitle = header.na(id(n,0),col(n,0))

                loc2 = where((stime_b(0,n,*) gt 0) and 			$
	                     (sdata_b(0,n,*) ne -1.0e32), count2)

                if (ind eq 0) then begin

		  if count-1 eq 0 then ystyle = 0 else ystyle = -8

                  plot, stime_b(0,n,loc2)-basetime, 			$
		        sdata_b(0,n,loc2),				$
                        xtickname = xtickname_temp, 			$
                        xtickv = xtickvalue, xticks = xtickn, 		$
                        xminor = xminor, xtitle=xtitle_temp,		$
	                xrange = [srtime, ertime],			$
                        xstyle = 1, pos = pos, /noerase,		$
	                ytitle = ytitle,				$
			ystyle = ystyle

		endif else begin

                  plot, stime_b(0,n,loc2)-basetime, 			$
		        sdata_b(0,n,loc2),				$
                        xstyle = 5, pos = pos, /noerase,		$
			ystyle = 5,					$
		        linestyle = 2
		  axis, yaxis=1,ytitle = ytitle

                endelse

              endfor

	      plot_num = plot_num + 1

            endif

          endif

          alt = alts(i)

          pressure = float(alt)*(log(yrange(1)) - log(yrange(0)))/float(ny) + $
               log(yrange(0))
          pressure = 'Pressure : '+string(10.0^pressure,format = '(e8.2)')

          get_position, ppp,space, sizes, plot_num, pos, /rect, 	$
		    xmargin = xmargin, ymargin=ymargin

          pos(0) = pos(0) - xmargin + space*3.0
          pos(2) = pos(2) - xmargin + space*3.0

          if plot_num eq ppp-1 then begin
            xtickname_temp = xtickname
            xtitle_temp = xtitle
          endif else begin
            xtickname_temp = xtickname2
            xtitle_temp = xtitle2
          endelse

          ytitle = header.na(id(var,1),col(var,1))

          if take_sd(var) then ytitle = strcompress(ytitle) + ' (Std. Dev.)'
          if take_mn(var) then ytitle = strcompress(ytitle) + ' (Mean)'

	  loc = where(binned_save(j,var,*,alt) ne 0.0,count)

          if count gt 0 then begin

	    if logged(var) eq 0 then 					$
              plot_io, time(loc), binned_save_mn(j,var,loc,alt), 	$
                xtickname = xtickname_temp, 				$
                xtickv = xtickvalue, xticks = xtickn, 			$
                xminor = xminor, xtitle=xtitle_temp,			$
	        xrange = [srtime, ertime],		        	$
                xstyle = 1, pos = pos, /noerase,			$
	        ytitle = ytitle						$
	    else							$
              plot, time(loc), binned_save_mn(j,var,loc,alt), 		$
                xtickname = xtickname_temp, 				$
                xtickv = xtickvalue, xticks = xtickn, 			$
                xminor = xminor, xtitle=xtitle_temp,			$
	        xrange = [srtime, ertime],		        	$
                xstyle = 1, pos = pos, /noerase,			$
	        ytitle = ytitle

;	    oplot, time(loc), binned_save(j,var,loc,alt), linestyle = 1

	    xyouts, pos(0)+0.01,pos(1)+0.005, pressure, charsize=0.8, /norm

          endif

	  plot_num = (plot_num + 1) mod ppp

        endfor

      endfor

    endfor

  endif

endif

if (!d.name eq 'PS') then begin

  device, /close
  set_plot, 'X'

endif

end
