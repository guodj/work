PRO fft_pro, single, t, dt, numofvar, numofpts, titlenam, unitnam,	      $
	data, locator, loc_time, title, unit, basetime

  maxi = 0

  for i=0,numofvar do begin

    loc = where(dt(i,*) ne 0.0, count)

    if count gt 0 then begin

        mini = min(dt(i,loc))
        steps = fix(max(t(i,*))/mini)+1

        if steps gt maxi then maxi=steps

    endif

  endfor

  for i=0, numofvar do begin

    loc2 = where(dt(i,*) gt 3.0*maxi, count)

    if count gt 0 then numofpts(i) = -1

  endfor

  if n_elements(maxi) eq 0 then maxi = max(numofpts)

  loc = where(numofpts gt -1, count)

  if count gt 0 then begin

    data_dum = fltarr(numofvar+1, maxi)

    for i=0,count-1 do begin

      ii = loc(i)

      loc2 = where(dt(ii,*) ne 0.0)

      mini = min(dt(ii,loc2))

      misdata = where(dt(ii,*) gt mini, count3)
      
      if count3 eq 0 then 						      $
	data_dum(i,0:numofpts(ii)) = single(ii,0:numofpts(ii))    $
      else begin

	npts = 0
        spt = 0

        for j=0,count3 do begin

	  if j gt 0 then spt = misdata(j-1)+1
	  if j ne count3 then						      $
	    ept = misdata(j)						      $
	  else								      $
	    ept = numofpts(ii)

	  dspt = npts
	  dept = npts + ept - spt

	  data_dum(ii,dspt:dept) = single(ii,spt:ept)

	  if j ne count3 then begin

	    npts = dept + 1
	    ept = dt(misdata(j))/mini

	    for k = 1, ept-1 do begin

	      data_dum(ii,npts) = (float(k)/float(ept))*		      $
		(single(ii,misdata(j))-single(ii,misdata(j)+1)) + 	      $
	        single(ii,misdata(j)+1)

	      npts = npts + 1

	    endfor

	  endif

	endfor

        numofpts(ii) = npts-1

      endelse

    endfor

  endif

  loc = where(numofpts ne max(numofpts),count)

  if count gt 0 then numofpts(loc) = -1

  loc = where(numofpts ne -1, count)

  if count gt 0 then begin

    data = fltarr(numofvar+3, max(numofpts)+1)

    for i=0,count-1 do data(loc(i)+2,*) = 				      $
      2.0*abs(fft(data_dum(loc(i),*),1))^2.0

    loc2 = where(dt ne 0.0)

    mini = min(dt(loc2))

    nt = float(max(numofpts)+1)*mini

    print, 'NT = ',nt, mini, max(numofpts)+1

    for i=0,max(numofpts)/2 do begin

      data(0,i) = float(i+1)/nt
      data(1,i) = 1.0/data(0,i)

    endfor

    data = data(*,0:max(numofpts)/2)

    numofpts(loc) = max(numofpts)/2

    numofpts = [max(numofpts),max(numofpts),numofpts]
    locator = intarr(max(numofpts)+1,2,numofvar+3)-1
    loc_time = dblarr(max(numofpts)+1)
    numofvar = numofvar + 2

    t_stamp = min(basetime) + double((min(t)+max(t))/2.0)

    for i=0,max(numofpts) do begin

      locator(i,0,loc+2) = loc+2
      locator(i,1,loc+2) = i
      locator(i,0,0) = 0
      locator(i,1,0) = i
      locator(i,0,1) = 1
      locator(i,1,1) = i
      loc_time(i) = t_stamp

    endfor

    title = ['Frequency','Period', titlenam]
    unit  = ['Hertz','Seconds',strarr(numofvar+1)+'dB']      

  endif

  return

end

