function amean, data, missing, log=log

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

    if n_e gt 1 then mn = mean(d(loc)) else mn = d(loc(0))

    if n_elements(log) gt 0 then mn = 10.0^mn

  endif else mn = missing

  return, mn

end

