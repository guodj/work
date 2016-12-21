
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = findfile('-t b*')
iFile = 0
while (strpos(initial_guess(iFile),"data") gt 0 or $
       strpos(initial_guess(iFile),"ps") gt 0 or $
       strpos(initial_guess(iFile),"sum") gt 0) do iFile = iFile + 1
initial_guess = initial_guess(iFile)
test = findfile(initial_guess)
if strlen(test(0)) eq 0 then initial_guess=initial_guess+'*'

filelist = ask('AMIE binary file name',initial_guess)

amiefiles = findfile(filelist)
nfiles = n_elements(amiefiles)

if (nfiles ge 1 and strlen(amiefiles(0)) gt 0) then begin

  for nf = 0, nfiles-1 do begin

    amie_file = amiefiles(nf)

    print, "reading ",amie_file

    ppp = 7

    read_amie_binary, amie_file, data, lats, mlts, time, fields, $
                  imf, ae, dst, hp, cpcp, version

    p = strpos(amie_file,"save")
    if (p gt 0) then amie_file = strmid(amie_file,0,p-1)

    psfile = amie_file+'_flat.ps'

    dlat = (lats(0) - lats(1)) * !dtor
    dlon = (mlts(1) - mlts(0)) * !pi / 12.0

    read_amie_data, amie_file+"_data", $
      RawMagData, AMIEMagData, MagFlag, MinFlag, $
      MagTime, MagList, nMags, nMagsTotal

    ;--------------------------------------------------------------
    ; figure out grid:
    ;--------------------------------------------------------------

    nlats = n_elements(lats)
    nmlts = n_elements(mlts)
    ntimes = n_elements(time)

    geoarea = fltarr(nmlts,nlats)

    re = 6372.0*1000.0

    for i=0,nlats-1 do for j=0,nmlts-1 do begin
      geoarea(j,i) = cos(lats(i)*!dtor) * dlat * dlon * re * re
    endfor

    nfields = n_elements(fields)

    type_1 = 6
    type_2 = 12

    ;--------------------------------------------------------------
    ; Get start time and end time, with defaults as the file times
    ;--------------------------------------------------------------

    c_r_to_a, itime_start, min(time)
    c_r_to_a, itime_end, max(time)

    c_a_to_s, itime_start, stime_start
    c_a_to_s, itime_end, stime_end

    start_time = stime_start

    c_s_to_a, itime_start, start_time
    c_a_to_r, itime_start, rtime_start

    end_time   = stime_end

    ;--------------------------------------------------------------
    ; If the user entered a short string, assume it is just a time
    ; and add the date on the front
    ;--------------------------------------------------------------

    if (strlen(end_time) lt 9) then begin
      if (strpos(end_time,'+') lt 0) then					$
    	end_time = strmid(start_time,0,9)+' '+end_time			$
      else begin
    	rtime_end = rtime_start+double(strmid(end_time,1,strlen(end_time)-1))*60.0
    	c_r_to_a, itime_end, rtime_end
    	c_a_to_s, itime_end, end_time
      endelse
    endif

    ;--------------------------------------------------------------
    ; Now figure out where in the file these things are, with the
    ; default to give the user everything
    ;--------------------------------------------------------------

    c_s_to_a, itime_end, end_time
    c_a_to_r, itime_end, rtime_end

    n_start = where(time ge rtime_start)
    if n_start(0) ge 0 then n_start = n_start(0) else n_start = 0

    n_end = where(time ge rtime_end)
    if n_end(0) ge 0 then n_end = n_end(0) else n_end = n_elements(time)-1

    ;--------------------------------------------------------------
    ; Put the contour data into data_1 array and get field name
    ;--------------------------------------------------------------

    data_1 = reform(data(*,type_1,*,*))
    data_2 = reform(data(*,type_2,*,*))
    field_1 = strcompress(fields(type_1))
    field_2 = strcompress(fields(type_2))

    ;--------------------------------------------------------------
    ; figure out contour levels:
    ;--------------------------------------------------------------

    maxi_array = fltarr(n_end-n_start+1)
    mini_array = fltarr(n_end-n_start+1)

    for i = n_start, n_end do begin
      mini_array(i-n_start) = min(data_1(i,*,*))
      maxi_array(i-n_start) = max(data_1(i,*,*))
    endfor


    mini_1  = min(data_1(n_start:n_end,*,*))
    maxi_1  = max(data_1(n_start:n_end,*,*))

    count = 0

    onepercent = n_elements(data_1(n_start:n_end,*,*)) * 0.0001

    while count lt onepercent do begin

      mini_1 = mini_1 * 0.95
      maxi_1 = maxi_1 * 0.95

      loc = where(data_1(n_start:n_end,*,*) le mini_1 or $
    		  data_1(n_start:n_end,*,*) ge maxi_1, count)

    endwhile

    if (mini_1 ge 0.0) then begin
      mini_1  = 0.0
      range_1 = maxi_1
    endif else range_1 = max([abs(mini_1),maxi_1])*2.0
    dc      = 10.0^fix(alog10(range_1/20.0))
    factor  = 0.0
    while (range_1 gt dc*20.0*factor) do factor=factor+0.05
    dc = factor*dc
    ;dc = float(ask('contour level for '+field_1,strcompress(string(dc))))
    if (mini_1) eq 0.0 then 						$
      levels_1 = findgen(21)*dc						$
    else levels_1 = (findgen(21) - 10.0)*dc

    ;--------------------------------------------------------------
    ; Set up device
    ;--------------------------------------------------------------

    setdevice, psfile, 'p', 5, 0.95

    plotdumb

    ;--------------------------------------------------------------
    ; Read color table. Blue to white to red for -/+ data,
    ; red to white for + only data.
    ;--------------------------------------------------------------

    makect, 'mid'
    ncolors = 255
    clevels = (ncolors-1)*findgen(21)/20.0 + 1

    ;--------------------------------------------------------------
    ; Set up plot sizes for the output graphs
    ;--------------------------------------------------------------

    space = 0.01
    pos_space, ppp, space, sizes, ny = ppp

    ;--------------------------------------------------------------
    ; Determine the character size in normalized coordinates
    ;--------------------------------------------------------------

    dy  = float(!d.y_ch_size)/float(!d.y_size)

    ;--------------------------------------------------------------
    ; If we have a whole bunch of plots, scale the size down.
    ;--------------------------------------------------------------

    if (ppp gt 12) then ch_size = 0.65 else ch_size = 1.0

    ntimes = n_end - n_start + 1

    image = fltarr(ntimes, nmlts)
    latimage = fltarr(ntimes, nmlts)
    area = fltarr(ntimes)
    cpcp = fltarr(ntimes)
    cpcpm = fltarr(ntimes)
    efield = fltarr(ntimes)
    efe = reform(data(*,7,*,*))
    efn = reform(data(*,8,*,*))
    eft = sqrt(efn^2 + efe^2)
    sigmap = reform(data(*,1,*,*))
    JouleHeating = SigmaP * Eft*Eft

    dl = (nmlts-1)/4
    ef_dayside = reform(eft(*,*,0:4))

    dl = 24.0 / (nmlts-1) * 15.0

    jh = fltarr(ntimes)

    for n = n_start, n_end do begin
      nn = n-n_start
      area(nn) = 0.0
      efield(nn) = max(ef_dayside(n,*,*))

      cpcp(nn) = max(data(n,0,*,*)) - min(data(n,0,*,*))
      cpcpm(nn) = max(data(n,17,*,*)) - min(data(n,17,*,*))

      hp(nn,1) = 0.0
      jh(nn) = 0.0

      for j=0,nmlts-1 do begin

    	; compute HPI

    	hp(nn,1) = hp(nn,1) +total(data_1(n,j,*)*geoarea(j,*)/1000)
    	jh(nn)   = jh(nn) + total(JouleHeating(n,j,*)*geoarea(j,*))

    	maxi = max(abs(data_1(n,j,*)))
    	mini = min(data_1(n,j,*))
    	if (abs(mini) eq maxi) then image(nn,j) = mini else image(nn,j) = maxi
    	loc = where(data_1(n,j,*) eq image(nn,j))
        if (loc(0) eq -1) then loc(0) = nlats/2
    	latimage(nn,j) = lats(loc(0))

    	ilat = loc(0)-1
    	done = 0
    	while (not done) do begin
    	  if (ilat le 0) then begin
    	    done = 1
    	  endif else begin
    	    if (abs(data_1(n,j,ilat)) gt 0.25*abs(data_1(n,j,loc(0)))) then begin
    	      ilat = ilat - 1
    	    endif else done = 1
    	  endelse
    	endwhile

    	if (ilat ge 0) then l = lats(ilat) else l = 90.0

    	area(nn) = area(nn) + (90.0 - l)*111.0 * dl * 111.0 * sin((90.0 - l)*!dtor)
      endfor
    endfor

    outputfile = amie_file+'_sum'

    OPENW,2,outputfile
    printf,2,' yyyy mm dd hh mm ss'+ $
    	     '      CPCP'+ $
    	     '        Ae'+ $
    	     '        Au'+ $
    	     '        Al'+ $
    	     '        Bx'+ $
    	     '        By'+ $
    	     '        Bz'+ $
    	     '       Vsw'+ $
    	     '  Comp Dst'+ $
    	     '  Comp HPI'+ $
    	     '  Comp SJH'+ $
    	     '   A (m^2)'+ $
             '      Emax'+ $
             '      nMag'

    for i =0,n_elements(time)-1 do begin
     c_r_to_a, itime, time(i)
     printf,2,format='(i5,5i3,13e10.2,i4)', itime, cpcp(i),$
    	      ae(i,0:2), imf(i,*), dst(i,0), hp(i,1),  $
              jh(i)/1000.0, area(i)*1.0e6, $
              efield(i), nMags(i)
    endfor

    CLOSE,2

    time_axis, rtime_start, rtime_end, btr, etr, $
    	       xtickname, xtitle, xtickv, xminor, xtickn

    get_position, ppp, space, sizes, 0, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    image = shift(image,0,-(nmlts-1)/2)
    image = alog10(image)

    latimage = shift(latimage,0,-(nmlts-1)/2)

    mini_1 = min(image)
    maxi_1 = max(image)
    range_1 = maxi_1-mini_1
    image = (image - mini_1)/range_1 * float(ncolors-1)

    mm_latimage = mm(latimage)

    latimage = (latimage - min(latimage)) / $
    	       (max(latimage)-min(latimage)) * float(ncolors-1)

;    tv,image, pos(0)+0.0005,pos(1)+0.0005,xsize=xsize,ysize=ysize,/norm
;
;    plot, [btr,etr], [0.0,1.0], /nodata, $
;    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
;    	  xtickname = strarr(30)+' ',	$
;    	  xtickv = xtickv,			$
;    	  xminor = xminor,			$
;    	  xticks = xtickn, xrange = [btr,etr],  $
;    	  ytickv = [0.0, 0.25, 0.5, 0.75, 1.0], $
;    	  ytickname = ['12','18','24','06','12'], yticks = 4, yminor = 6, $
;    	  ytitle = 'MLT'
;    ctpos = [pos(2)+0.01, pos(1), pos(2)+0.04, pos(3)]
;    ctmm = [mini_1, maxi_1]
;    cttitle = 'log('+field_1+')'
;    plotct, ncolors, ctpos, ctmm, cttitle, /right


    xyouts, 1.0, 1.01, "AMIE Version "+string(Version,format="(f5.2)"), $
      alignment = 1.0, /norm

    !p.thick = 3


    get_position, ppp, space, sizes, 0, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    plot, time(n_start:n_end)-rtime_start, alog10(hp(*,1)), $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ',	$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr], $
      ytitle = 'log(HPI)'

    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'log(Watts)' , $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 1, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    plot, time(n_start:n_end)-rtime_start, alog10(jh/1000.0), $
    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
    	  xtickname = strarr(30)+' ',	$
    	  xtickv = xtickv,			$
    	  xminor = xminor,			$
    	  xticks = xtickn, xrange = [btr,etr], $
      ytitle = 'log(SJH)'

    xyouts, pos(2)+0.04, (pos(1)+pos(3))/2.0, 'Joule Heating' , $
    	    orientation = 270, alignment = 0.5, /norm
    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'log(Watts)' , $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 2, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    plot, time(n_start:n_end)-rtime_start, area/1.0e6, $
    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
    	  xtickname = strarr(30)+' ',	$
    	  xtickv = xtickv,			$
    	  xminor = xminor,			$
    	  xticks = xtickn, xrange = [btr,etr], $
      ytitle = 'Area'

    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, '(m!E2!N/1.0e12)', $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 3, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    mini = min([floor(min(ae(n_start:n_end,2))/100.0)*100.0, -500.0])
    maxi = max([fix(max(ae(n_start:n_end,2))/10.0+1.0)*10.0, 0.0])

    plot, time(n_start:n_end)-rtime_start, ae(n_start:n_end,2), $
    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
    	  xtickname = strarr(30)+' ',	$
    	  xtickv = xtickv,			$
    	  xminor = xminor,			$
    	  xticks = xtickn, xrange = [btr,etr], yrange=[mini,maxi], $
      ytitle = 'AL'

    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, '(nT)', $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 4, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    mini = min([floor(min(dst(n_start:n_end,0))/50.0)*50.0, -50.0])
    maxi = max([fix(max(dst(n_start:n_end,0))/10.0+1.0)*10.0, 0.0])

    plot, time(n_start:n_end)-rtime_start, dst(n_start:n_end,0), $
    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
    	  xtickname = strarr(30)+' ',	$
    	  xtickv = xtickv,			$
    	  xminor = xminor,			$
    	  xticks = xtickn, xrange = [btr,etr], yrange=[mini,maxi], $
      ytitle = 'Dst'

    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, '(nT)', $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 5, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    mini = 0.0
    maxi = max([fix(max(cpcp)/50.0+1.0)*50.0, 100.0])

    print, "mean cpcp : ", mean(cpcp)

    plot, time(n_start:n_end)-rtime_start, cpcp, $
    	  xstyle = 1, ystyle = 1, pos = pos, /noerase, $
    	  xtickname = strarr(30)+' ',	$
    	  xtickv = xtickv,			$
    	  xminor = xminor,			$
    	  xticks = xtickn, xrange = [btr,etr], yrange=[mini,maxi], $
      ytitle = 'CPCP'

    oplot, time(n_start:n_end)-rtime_start, cpcpm, linestyle = 2

    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, '(kV)', $
    	    orientation = 270, alignment = 0.5, /norm

    get_position, ppp, space, sizes, 6, pos, /rect
    pos(2) = pos(2) - 0.05
    pos(0) = pos(0) + 0.05
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001

    plot, MagTime(n_start:n_end)-rtime_start, nMags(n_start:n_end), $
	   xstyle = 1, ystyle = 1, pos = pos, /noerase, $
	   xtickname = xtickname,			$
	   xtitle = xtitle,			$
	   xtickv = xtickv,			$
	   xminor = xminor,			$
	   xticks = xtickn, xrange = [btr,etr], $
      ytitle = 'nMags', $
           yrange = [min(nMags(n_start:n_end))-1,max(nMags(n_start:n_end))+1]

    xyouts, pos(2)+0.04, (pos(1)+pos(3))/2.0, 'Number of' , $
    	    orientation = 270, alignment = 0.5, /norm
    xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'Magnetometers' , $
    	    orientation = 270, alignment = 0.5, /norm

;    plot, time(n_start:n_end)-rtime_start, imf(n_start:n_end,2), $
;	   xstyle = 1, ystyle = 1, pos = pos, /noerase, $
;	   xtickname = xtickname,			$
;	   xtitle = xtitle,			$
;	   xtickv = xtickv,			$
;	   xminor = xminor,			$
;	   xticks = xtickn, xrange = [btr,etr],  $
;	   yrange = mm(imf(n_start:n_end,1:2))
;
;    oplot, time(n_start:n_end)-rtime_start, imf(n_start:n_end,1), linestyle = 2

    oplot, [btr,etr], [0.0,0.0], linestyle = 1

    ;xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'E-Field (mV/m)', $
	     ;orientation = 270, alignment = 0.5, /norm

    closedevice

  endfor

endif

end





