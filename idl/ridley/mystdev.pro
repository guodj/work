function stdev, data, missing

  if n_elements(missing) eq 0 then missing = 2.0*max(data)

  loc = where(data ne missing, count)
  if count gt 1 then begin
    mn = mean(data(loc))
    diff = (data(loc)-mn)^2.0
    sum = 0.0
    for n=0,count-1 do sum = sum + diff(n)
    sum = (sum/float(count-1))^0.5
  endif else sum = missing
  return, sum

end

