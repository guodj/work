;
; procedure det
;
;  Determines the determininant of a square array recursively
;

function det, array

  if n_elements(array(*,0)) eq 2 then begin
    value = array(0,0)*array(1,1) - array(1,0)*array(0,1)
  endif else begin

    for i=0,n_elements(array(0,*))-1 do begin

      array2 = fltarr(n_elements(array(0,*))-1, n_elements(array(*,0))-1)

      for j=0,i-1 do 							$
	array2(j,0:n_elements(array(*,0))-2) =				$
	  array(j,1:n_elements(array(*,0))-1)

      for j=i+1,n_elements(array(0,*))-1 do 				$
	array2(j-1,0:n_elements(array(*,0))-2) =			$
	  array(j,1:n_elements(array(*,0))-1)

      value2 = 0.0

      value2 = det(array2)

      if i mod 2 eq 0 then fac = +1.0 else fac = -1.0

      value = value + fac*value2*array(i,0)

    endfor

  endelse

  return, value

end
