pro smooth, x,y, missing 

  nomiss = where(x ne missing, nocount)
  loc = where(x eq missing, count)
  if (count gt 0) and (nocount gt 0) then begin
    x(loc) = min(x(nomiss))
    y(loc) = 100*max(y(nomiss))
  endif

  nx = n_elements(x(*,0))
  ny = n_elements(x(0,*))

  for i=0,nx-1 do for j=0,ny-1 do begin

; compute weights from error terms

    if j gt 0 then begin
      if y(i,j-1) gt 0.0 then ry1 = 1.0/y(i,j-1) else 		$
         ry1 = 1.0e-4
    endif
    if y(i,j) gt 0.0 then ry2 = 1.0/y(i,j) else 		$
       ry2 = 1.0e-4
    if j lt ny-1 then begin
      if y(i,j+1) gt 0.0 then ry3 = 1.0/y(i,j+1) else 		$
       ry3 = 1.0e-4
    endif

    if i gt 0 then begin
      if y(i-1,j) gt 0.0 then rx1 = 1.0/y(i-1,j) else 		$
        rx1 = 1.0e-4
    endif
    if y(i,j) gt 0.0 then rx2 = 1.0/y(i,j) else 		$
      rx2 = 1.0e-4
    if i lt nx-1 then begin
      if y(i+1,j) gt 0.0 then rx3 = 1.0/y(i+1,j) else 		$
        rx3 = 1.0e-4
    endif

    if rx2 gt 1.0e-4 then begin

; Center section

      if (j gt 0) and (j lt ny-1) and				$
         (i gt 0) and (i lt nx-1) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx1*x(i-1,j) + rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry1+ry2+ry3+rx1+rx2+rx3)

      endif

; front side

      if (j eq 0) and (j ne ny-1) and				$
         (i gt 0) and (i lt nx-1) then begin

        x(i,j) = (ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx1*x(i-1,j) + rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry2+ry3+rx1+rx2+rx3)

      endif

; Back side

      if (j eq ny-1) and (j ne 0) and				$
         (i gt 0) and (i lt nx-1) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) +	$
	         rx1*x(i-1,j) + rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry1+ry2+rx1+rx2+rx3)

      endif

; Left side

      if (j gt 0) and (j lt ny-1) and				$
         (i eq 0) and (i ne nx-1) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry1+ry2+ry3+rx2+rx3)

      endif

; Right side

      if (j gt 0) and (j lt ny-1) and				$
         (i eq nx-1) and (i ne 0) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx1*x(i-1,j) + rx2*x(i,j))/	$
	         (ry1+ry2+ry3+rx1+rx2)

      endif

; front left corner

      if (j eq 0) and (j ne ny-1) and				$
         (i eq 0) and (i ne nx-1) then begin

        x(i,j) = (ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry2+ry3+rx2+rx3)

      endif

; back right corner

      if (j eq ny-1) and (j ne 0) and				$
         (i eq nx-1) and (i ne 0) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) +	$
	         rx2*x(i,j) + rx1*x(i-1,j))/	$
	         (ry1+ry2+rx2+rx1)

      endif

; back left corner

      if (j eq ny-1) and (j ne 0) and				$
         (i eq 0) and (i ne nx-1) then begin

        x(i,j) = (ry1*x(i,j-1) + ry2*x(i,j) +	$
	         rx2*x(i,j) + rx3*x(i+1,j))/	$
	         (ry1+ry2+rx2+rx3)

      endif

; front right corner

      if (j eq 0) and (j ne ny-1) and				$
         (i eq nx-1) and (i ne 0) then begin

        x(i,j) = (ry2*x(i,j) + ry3*x(i,j+1) +	$
	         rx1*x(i-1,j) + rx2*x(i,j))/	$
	         (ry2+ry3+rx1+rx2)

      endif

; if only 1d array :

      if (j eq 0) and (j eq ny-1) then begin

        if (i eq 0) and (i ne nx-1) then begin

	  x(i) = (rx2*x(i) + rx3*x(i+1) + ry2*x(i))/	$
		 (rx2 + rx3 + ry2)

	endif

        if (i ne 0) and (i ne nx-1) then begin

	  x(i) = (rx1*x(i-1) + rx2*x(i) + rx3*x(i+1) + ry2*x(i))/	$
		 (rx1 + rx2 + rx3 + ry2)

	endif

        if (i ne 0) and (i eq nx-1) then begin

	  x(i) = (rx1*x(i-1) + rx2*x(i) + ry2*x(i))/	$
		 (rx1 + rx2 + ry2)

	endif

      endif

    endif

  endfor

  return

end
