
PRO sort, array, order, arraymd

  n = n_elements(array) - 1

  sorted = dblarr(n+1)
  order = lonarr(n+1)

  for i=0L,n do begin

    find = where(array lt array(i), count)

    sorted(count) = array(i)
    order(count) = i

  endfor

  array = sorted

  if n_elements(arraymd) gt 0 then begin

    arraymd(0) = array(0)

    npts = 0L

    for i=1L,n do begin

      loc = where(arraymd eq array(i), count)

      if count eq 0 then begin

	npts = npts + 1
	arraymd(npts) = array(i)

      endif

    endfor

    arraymd = arraymd(0:npts)

  endif    

END
