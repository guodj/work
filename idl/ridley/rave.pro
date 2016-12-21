PRO rave, avearr, nadd

;	How many points are there in the array:

  npts = n_elements(avearr)-1

;	If there is enough data, go ahead with the filter

  if (npts gt nadd) and (nadd ge 3) then begin

;	We need another array, so we don't double count points. In other words,
;	filter the old array, while shoving the filtered values into a new
;	array.

    dumarr = fltarr(npts+1)

;	Since the first 4 and last 4 points are not filtered, just put them
;	into the new array

    left = (nadd-1)/2
    right = nadd - left - 1

    dumarr(0:left-1) = avearr(0:left-1)
    dumarr(npts-right-1:npts) = avearr(npts-right-1:npts)

;	Filter the rest

    for i = left, npts-right do begin

      sum = 0.0

      for j=-left,right do sum=sum+avearr(i+j)

      dumarr(i) = sum / nadd

    endfor

;	Put the new data beck into the original data array

    avearr = dumarr

  endif

  RETURN

END
