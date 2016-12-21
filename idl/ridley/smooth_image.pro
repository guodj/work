pro smooth_image, image, missing

  xn = n_elements(image(*,0))
  yn = n_elements(image(0,*))
  if yn eq 0 then yn = 1

  for i=0,xn-1 do for j=0,yn-1 do begin

    total = 0.0
    weight = 0.0

    if image(i,j) eq missing then begin
      main_i = 0.0
      main_weight = 0.0
    endif else begin
      main_i = image(i,j)
      main_weight = 1.0
    endelse

    total = total + main_i
    weight = weight + main_weight

    pt = i-1
    pta = 0.0
    pta_weight = 0.0

    while pt ge 0 and pta eq 0.0 do begin
      if image(pt,j) ne missing then begin
	pta_weight = 1.0/(float(i-pt)*2.0)
        pta = image(pt,j)*pta_weight
	if (image(pt,j) eq 0.0) then pta = 1.0e-6
      endif else pt = pt - 1
    endwhile

    total = total + pta
    weight = weight + pta_weight

    pt = i+1
    pta = 0.0
    pta_weight = 0.0

    while pt le xn-1 and pta eq 0.0 do begin

      if image(pt,j) ne missing then begin
	pta_weight = 1.0/(float(pt-i)*2.0)
        pta = image(pt,j)*pta_weight
	if (image(pt,j) eq 0.0) then pta = 1.0e-6
      endif else pt = pt + 1
    endwhile

    total = total + pta
    weight = weight + pta_weight

    pt = j-1
    pta = 0.0
    pta_weight = 0.0

    while pt ge 0 and pta eq 0.0 do begin
      if image(i,pt) ne missing then begin
	pta_weight = 1.0/(float(j-pt)*2.0)
        pta = image(i,pt)*pta_weight
	if (image(i,pt) eq 0.0) then pta = 1.0e-6
      endif else pt = pt - 1
    endwhile

    total = total + pta
    weight = weight + pta_weight

    pt = j+1
    pta = 0.0
    pta_weight = 0.0

    while pt le yn-1 and pta eq 0.0 do begin
      if image(i,pt) ne missing then begin
	pta_weight = 1.0/(float(pt-j)*2.0)
        pta = image(i,pt)*pta_weight
	if (image(i,pt) eq 0.0) then pta = 1.0e-6
      endif else pt = pt + 1
    endwhile

    total = total + pta
    weight = weight + pta_weight

    if weight gt 0.0 then image(i,j) = total/weight

  endfor

  return

end