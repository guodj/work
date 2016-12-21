pro sort_a, x, y

  nele = n_elements(x)
  nx = x
  ny = y
  maxi = max(x)*2.0

  for i=0,nele-1 do begin
    loc = where(nx eq min(nx))
    x(i) = nx(loc(0))
    y(i) = ny(loc(0))
    nx(loc(0)) = maxi
  endfor

  return

end
