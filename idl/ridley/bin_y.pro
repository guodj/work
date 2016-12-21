
pro bin_y, data, xaxis, nbin, output, smooth=smooth, range=range, 	$
	logbin=logbin, missing = missing, rotate = rotate

  if n_elements(range) lt 2 then begin
    mini = min(xaxis)
    maxi = max(xaxis)
  endif else begin
    mini = range(0)
    maxi = range(1)
  endelse

  if n_elements(missing) eq 0 then missing = 0.0

  if n_elements(smooth) gt 0 then nsm = smooth

  if n_elements(rotate) gt 0 then rotate = 1 else rotate = 0

  nalt = n_elements(data(0,*))
  nold = n_elements(xaxis)

  output = fltarr(nbin,nalt)
  npts = intarr(nbin,nalt)

  xnew = double(findgen(nbin))*(maxi-mini)/double(nbin) + mini

  if rotate then begin

    xaxis = xaxis mod 360.0

    loc = where(xaxis gt maxi, count)
    if count gt 0 then xaxis(loc) = xaxis(loc) - 360.0

    loc = where(xaxis lt mini, count)
    if count gt 0 then xaxis(loc) = xaxis(loc) + 360.0

  endif

  locrange = where(xaxis ge mini and xaxis le maxi, count)

  if count gt 0 then begin

    i = long(0)

    for i=long(0),long(count-1) do begin

      ii = locrange(i)

      d = abs(xaxis(locrange(i)) - xnew)
      loc = where(d eq min(d))

      if n_elements(logbin) gt 0 then begin
        for j=0,nalt-1 do begin
          if data(ii,j) gt 0.0 then begin
	    output(loc(0),j) = output(loc(0),j) + log(data(ii,j))
	    npts(loc(0),j) = npts(loc(0),j) + 1
	  endif
	endfor
      endif else begin
        for j=0,nalt-1 do begin
	  output(loc(0),j) = output(loc(0),j) + data(ii,j)
	  npts(loc(0),j) = npts(loc(0),j) + 1
	endfor
      endelse

    endfor

    loc = where(npts gt 0,count)
    if count gt 0 then begin
      if n_elements(logbin) gt 0 then 					$
        output(loc) = 10.0^(output(loc)/npts(loc))			$
      else								$
        output(loc) = output(loc)/npts(loc)
    endif

    loc = where(npts eq 0,count)
    if count gt 0 then output(loc) = missing

    if n_elements(smooth) gt 0 then begin

      for i=1,nbin-2 do for k=0,nalt-1 do begin

        if output(i,k) eq missing then begin

          j = i+1
          while ((output(j,k) eq missing) and 				$
		 (j lt i+nsm) and 					$
		 (j lt nbin-1)) do j=j+1

          if (output(j,k) ne missing) and 				$
	     (output(i-1,k) ne missing) then begin

            if n_elements(logbin) gt 0 then 				$
              output(i,k) = 10.0^((log(output(i-1,k)) + 		$
				   log(output(j,k))/float(j-i))/	$
				   (1.0+1.0/float(j-i)))		$
	    else							$
              output(i,k) = (output(i-1,k) + 			$
			     output(j,k)/float(j-i))/		$
			     (1.0+1.0/float(j-i))
	  endif

        endif

      endfor

    endif

  endif

  return

end