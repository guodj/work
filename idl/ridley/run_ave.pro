pro run_ave,data,npts,missing

  n = n_elements(data)

  for i=0,n-1 do begin

    bj = i-npts/2
    if bj lt 0 then bj = 0
    ej = bj + npts
    if ej gt n-1 then ej = n-1

    nave = 0.0
    ave = 0.0

    for j=bj,ej do begin

      if data(j) ne missing then begin
        ave = ave + data(j)  ; /float(abs(j-i)+1)
        nave = nave + 1.0    ; /float(abs(j-i)+1)
      endif

    endfor

    if nave gt 0.0 then data(i) = ave/nave

  endfor

  return

end
