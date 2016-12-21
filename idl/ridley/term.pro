function term, long, day, time

  tilt = 23.5/!radeg
  if n_elements(time) eq 0 then time = 12.0
  if n_elements(day) eq 0 then day = 365.0/2.0
  daytheta = (day-365.0/4.0)*2.0*!pi/365.0

  st = sin(tilt)
  ct = cos(tilt)
  sd = sin(daytheta)
  cd = cos(daytheta)

  longrad = long/!radeg
  sl = sin(longrad-daytheta+time/24.0*2.0*!pi+!pi/2)
  cl = cos(longrad-daytheta+time/24.0*2.0*!pi+!pi/2)

  phi = (atan(-sd*st/(cl*sd*ct + sl*cd)))*!radeg

  loc = where(phi lt 0.0, count)
  if count gt 0 then phi(loc) = -90-phi(loc)
  loc = where(phi ge 0.0, count)
  if count gt 0 then phi(loc) = 90-phi(loc)

  return, phi

end
