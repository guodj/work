
common ffinfo, header

stime = [95,01,10,00,00,00]
etime = [95,01,12,00,00,00]

if n_elements(data_vec) gt 0 then que = ask("new data (y/n) ?","y")	$
else que = 'y'

if que eq 'y' then begin

  print, 'Reading data...'

  filename = 'imf_95_test'

  col_vec = [41,42,43,44]

  read_flat_vector, stime, etime, col_vec, time, data_vec, nrows,	$
	filename = filename

  col_scal = [0,1,2,3,13,14,15,16,40,29,30,31,26,27,28]

  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,	$
	filename = filename

  nrows = nrows(0)

endif

plotdumb

ppp = 5
space = 0.03

pos_space,ppp,space,sizes, nx = 1

for col = 0, 4 do begin

  get_position, ppp, space, sizes, col, pos, xmargin = 0.075, /rect

  delay = fltarr(nrows)
  cross = fltarr(nrows)
  magdel = fltarr(nrows)
  magdel(*) = data_scal(10,*)
  data_cross = fltarr(nrows,61)

  if col lt 4 then data_cross(*,*) = data_vec(col,*,*)

  for i=0,nrows-1 do begin

    if data_vec(0,i,0) gt -1.0e32 then begin

      if col eq 4 then begin

;	for j=0,60 do begin

;	  maxi = 0.0

;	  for k=0,2 do begin

;	    temp = max(data_vec(k,i,*)) - mean(data_vec(k,i,*))
;	    if temp gt 0.5 and max(data_vec(k,i,*)) gt maxi then begin
;	      sav = k
;	      maxi = max(data_vec(k,i,*))
;	    endif

;          endfor

;	  data_cross(i,*) = data_vec(sav,i,*)

;	  data_cross(i,j) = max(data_vec(0:3,i,j))

;        data_cross(i,*) = ((1.0+data_vec(0,i,*))^2.0 + 		$
;			   (1.0+data_vec(1,i,*))^2.0 + 			$
;			   (1.0+data_vec(2,i,*))^2.0 + 			$
;			   (1.0+data_vec(3,i,*))^2.0)^0.5/2.0

        data_cross(i,*) = (data_vec(0,i,*) +	 		$
			   data_vec(1,i,*) + 			$
			   data_vec(2,i,*))/3.0

;        endfor

      endif

      cross(i) = max(data_cross(i,*))
      loc = where(data_cross(i,*) eq cross(i))
      diff = (float(loc(0)) - 30.0)*2.0
      delay(i) = data_scal(8,i) + diff

    endif else delay(i) = 10000.0

  endfor

  loc = where(magdel eq -1.0e32, count)
  if count gt 0 then magdel(loc) = 10000.0

  loc = where(cross lt 0.90,count)

  if count gt 0 then begin
    delay(loc) = 10000.0
    magdel(loc) = 10000.0
  endif

  diff = abs(delay - magdel)
  if count gt 0 then diff(loc) = 10000.0
  
  loc = where(diff lt 1000.0,count)
  if count gt 0 then print, tostr(col),'. Mean difference : ',mean(diff(loc))

  plot, delay, max_value = 1000.0, /noerase, pos = pos
  oplot, magdel, max_value = 1000.0, linestyle = 2
  oplot, diff, max_value = 1000.0, linestyle = 1

endfor

end
