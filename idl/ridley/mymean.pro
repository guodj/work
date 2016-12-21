function mean, data, missing, log=log

  d = data

  if n_elements(missing) eq 0 then missing = 2.0*max(data)

  if n_elements(log) gt 0 then begin
    loc = where(d eq missing or d eq 0.0,count)
    if count gt 0 then d(loc) = missing
    loc = where(d ne missing,count)
    if count gt 0 then d(loc) = alog10(d(loc))
  endif

  loc = where(d ne missing,n_e)

  if n_e gt 0 then begin

    n = 0
    mean = 0.0

    for i=0,n_e-1 do mean = mean + d(loc(i))

    mean=mean/float(n_e)

    if n_elements(log) gt 0 then mean = 10.0^mean

  endif else mean = missing

  return, mean

end

