;
; procedure expand_array
;
;  expands a 2 dimensional array by a factor of SIZE
;

pro expand_array, array, fsize

  nx = n_elements(array(*,0))
  ny = n_elements(array(0,*))

  if n_elements(fsize) gt 1 then begin
    sizex = fsize(0)
    sizey = fsize(1)
  endif else begin
    sizex = fsize
    sizey = fsize
  endelse

  nx2 = nx*sizex - sizex + 1
  ny2 = ny*sizey - sizey + 1

  if ny gt 0 then tarray = fltarr(nx2,ny2) 		$
  else tarray = fltarr(nx2)

  if ny gt 0 then begin

    for i = 0, nx2 - 1, sizex do begin

      ii = i/sizex

      for j = 0, ny-2 do for k=0,sizey-1 do begin

        tarray(i,j*sizey+k) = 						$
	  float(k)*(array(ii,j+1)-array(ii,j))/float(sizey) + array(ii,j)

      endfor

      tarray(i,ny2-1) = array(ii, ny-1)

    endfor

  endif

  if ny eq 0 then ny2 = 1

  for j = 0, ny2 - 1 do begin

    for i = 0, nx-2 do for k=0,sizex-1 do begin

      tarray(i*sizex+k,j) = 						$
	float(k)*(tarray((i+1)*sizex,j)-tarray(i*sizex,j))/float(sizex) + $
	tarray(i*sizex,j)

    endfor

  endfor

  array = tarray

  return

end    

