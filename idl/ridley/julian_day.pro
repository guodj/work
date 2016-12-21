
function julian_day, time
  jd = fix(time*0.0)
  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  for i=0,n_elements(time)-1 do begin
      if time(i) gt 0 then begin
          c_r_to_a, itime, time(i)
          year = itime(0)
          mon = itime(1)
          day = itime(2)
          if year mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
          julianday = 0
          for i=0,mon-2 do julianday = julianday + dayofmon(i)
          jd(i) = julianday + day
      endif
  endfor
  return, jd
end
