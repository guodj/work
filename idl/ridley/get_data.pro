
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Procedure Get_Data
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO Get_Data, single, t, dt, basetime, numofvar, numofpts, plottype, 	$
	      times, datafile, grafvar, same

  begtime = fltarr(8)
  endtime = fltarr(8)
  datadum = fltarr(30000)
  timedum = fltarr(30000)
  basedum = double(0.0)
  numdum = 0.0

  testgra  = where(grafvar(0:numofvar,0) ne 0,cg)

  if (cg eq 0) and (same eq 1) then begin

    k = 0

    filelen = float(strlen(datafile(grafvar(k,0))))

    filedum = datafile(grafvar(k,0)) + $
	      '                                                       ' + $
	      '                                                       ' + $
	      '                                                       '
    filedum = strmid(filedum,0,100)

    begtime(0) = float(times(k,0,0)) + 1900.0
    begtime(1) = 0.0
    begtime(2) = float(times(k,0,1))
    begtime(3) = float(times(k,0,2))
    begtime(4) = float(times(k,0,3))
    begtime(5) = float(times(k,0,4))
    begtime(6) = float(times(k,0,5))
    begtime(7) = 0.0
    endtime(0) = float(times(k,1,0)) + 1900.0
    endtime(1) = 0.0
    endtime(2) = float(times(k,1,1))
    endtime(3) = float(times(k,1,2))
    endtime(4) = float(times(k,1,3))
    endtime(5) = float(times(k,1,4))
    endtime(6) = float(times(k,1,5))
    endtime(7) = 0.0

    bty=string(times(k,0,0))
    btm=string(times(k,0,1))
    btd=string(times(k,0,2))
    bth=string(times(k,0,3))
    if times(k,0,4) ge 10 then btmi=string(times(k,0,4)) else btmi='0'+string(times(k,0,4))
    if times(k,0,5) ge 10 then bts=string(times(k,0,5)) else bts='0'+string(times(k,0,5))
    ety=string(times(k,1,0))
    etm=string(times(k,1,1))
    etd=string(times(k,1,2))
    eth=string(times(k,1,3))
    if times(k,1,4) ge 10 then etmi=string(times(k,1,4)) else etmi='0'+string(times(k,1,4))
    if times(k,1,5) ge 10 then ets=string(times(k,1,5)) else ets='0'+string(times(k,1,5))

    print, 'Getting Data for interval:'
    print, 								    $
     strcompress(btm+'-'+btd+'-'+bty+','+bth+':'+btmi+':'+bts,/remove_all), $
     '  to  ',								    $
     strcompress(etm+'-'+etd+'-'+ety+','+eth+':'+etmi+':'+ets,/remove_all)

    dumnum = 200.0
    datadum = fltarr(long(numofvar+1)*long(numdum+1))
    timedum = fltarr(long(numofvar+1)*long(numdum+1))
    colnum = float(grafvar(0:numofvar,1)+1)

    z = call_external('readmulti_exe','read_flat', $
        filedum, begtime, endtime, datadum, timedum, basedum, $
        numdum, colnum, dumnum, filelen)

    datadum = fltarr(long(numofvar+1)*long(numdum+1))
    timedum = fltarr(long(numofvar+1)*long(numdum+1))
    t = fltarr(numofvar+1,numdum+1)
    dt = fltarr(numofvar+1,numdum+1)
    single = fltarr(numofvar+1,numdum+1)

    colnum = float(grafvar(0:numofvar,1)+1)

    dumnum = float(numofvar)

    if numdum gt 0 then begin	

      z = call_external('readmulti_exe','read_flat', $
          filedum, begtime, endtime, datadum, timedum, basedum, numdum, colnum, dumnum,filelen)

      for k=0, numofvar do begin

        basetime(k) = basedum

        for i=0, numdum do begin

          j = long(k) * long(numdum+1) + long(i)

          if (datadum(j) ne -1.0e32) or 				      $
	     ((plottype ge 8) and (plottype ne 20)) or 			      $
	     (plottype eq 5) then begin

	    if (datadum(j) ne -1.0e32) then 				      $
	      single(k,numofpts(k)) = datadum(j) 			      $
	    else 							      $
	      single(k,numofpts(k)) = -1.0e10
	    t(k,numofpts(k)) = timedum(j)
            numofpts(k) = numofpts(k) + 1

          endif

        endfor

        numofpts(k) = numofpts(k) - 1

        if numofpts(k) gt 1 then dt(k,0:numofpts(k)-2) = t(k,1:numofpts(k)-1) - t(k,0:numofpts(k)-2)

      endfor

    endif else begin

      print, 'No Data in time period selected!!!!'
      for k=0,numofvar do numofpts(k)=-1

    endelse

  endif else begin

    oldnum = 0

    for k=0, numofvar do begin

      filedum = datafile(grafvar(k,0)) + $
      	        '                                                       ' + $
	        '                                                       ' + $
	        '                                                       '
      filedum = strmid(filedum,0,100)

      begtime(0) = float(times(k,0,0)) + 1900.0
      begtime(1) = 0.0
      begtime(2) = float(times(k,0,1))
      begtime(3) = float(times(k,0,2))
      begtime(4) = float(times(k,0,3))
      begtime(5) = float(times(k,0,4))
      begtime(6) = float(times(k,0,5))
      begtime(7) = 0.0
      endtime(0) = float(times(k,1,0)) + 1900.0
      endtime(1) = 0.0
      endtime(2) = float(times(k,1,1))
      endtime(3) = float(times(k,1,2))
      endtime(4) = float(times(k,1,3))
      endtime(5) = float(times(k,1,4))
      endtime(6) = float(times(k,1,5))
      endtime(7) = 0.0

      bty=string(times(k,0,0))
      btm=string(times(k,0,1))
      btd=string(times(k,0,2))
      bth=string(times(k,0,3))
      if times(k,0,4) ge 10 then btmi=string(times(k,0,4)) else btmi='0'+string(times(k,0,4))
      if times(k,0,5) ge 10 then bts=string(times(k,0,5)) else bts='0'+string(times(k,0,5))
      ety=string(times(k,1,0))
      etm=string(times(k,1,1))
      etd=string(times(k,1,2))
      eth=string(times(k,1,3))
      if times(k,1,4) ge 10 then etmi=string(times(k,1,4)) else etmi='0'+string(times(k,1,4))
      if times(k,1,5) ge 10 then ets=string(times(k,1,5)) else ets='0'+string(times(k,1,5))

      colnum = float(grafvar(k,1)+1)

      z = call_external('readflat_exe','read_flat', $
          filedum, begtime, endtime, datadum, timedum, basedum, numdum, colnum)

      basetime(k) = basedum

      inumdum=fix(numdum)

      if k eq 0 then begin
	t = fltarr(numofvar+1,inumdum+1)
	dt = fltarr(numofvar+1,inumdum+1)
        single = fltarr(numofvar+1,inumdum+1)
	oldnum = inumdum
      endif

      if oldnum lt inumdum then begin
	tn=fltarr(numofvar+1,inumdum+1)
	dtn=fltarr(numofvar+1,inumdum+1)
	singlen=fltarr(numofvar+1,inumdum+1)
        tn(*,0:oldnum)=t(*,0:oldnum)
        dtn(*,0:oldnum)=dt(*,0:oldnum)
        singlen(*,0:oldnum)=single(*,0:oldnum)
	t=tn
	dt=dtn
	single=singlen
	tn=0
	dtn=0
	singlen=0
	oldnum=inumdum
      endif

      for i=0, inumdum do begin

        if (datadum(i) ne -1.0e32) or 					$
	   ((plottype ge 8) and						$ 
	    (plottype ne 20)) or					$
	   (plottype eq 5) then begin

	  if (datadum(i) ne -1.0e32) then 				      $
	    single(k,numofpts(k)) = datadum(i) 				      $
	  else 								      $
	    single(k,numofpts(k)) = -1.0e10
	  t(k,numofpts(k)) = timedum(i)
          numofpts(k) = numofpts(k) + 1

        endif

      endfor

      numofpts(k) = numofpts(k) - 1

      if numofpts(k) gt 1 then						      $
	dt(k,0:numofpts(k)-2) = t(k,1:numofpts(k)-1) - t(k,0:numofpts(k)-2)

    endfor

  endelse 

  if max(numofpts) gt 0 then begin
    single=single(*,0:max(numofpts))
    t=t(*,0:max(numofpts))
    dt=dt(*,0:max(numofpts)-1)
  endif else begin
    single=fltarr(numofvar+1,2)
    t=fltarr(numofvar+1,2)
    dt=fltarr(numofvar+1,2)
  endelse

  RETURN

END

