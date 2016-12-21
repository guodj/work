
pro smooth_1d, time, data, dt, outdata

  nTimes = n_elements(time)
  outData = fltarr(nTimes)

  for iTime = 0, nTimes-1 do begin
      loc = where(time ge time(iTime)-dt/2.0 and time lt time(iTime)+dt/2.0)
      OutData(iTime) = mean(data(loc))
  endfor

  return

end
