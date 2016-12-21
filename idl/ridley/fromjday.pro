function fromjday, year, jday
  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  if year mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
  mon = 1
  dum = jday
  while (dum gt dayofmon(mon-1)) do begin
    dum = dum - dayofmon(mon-1)
    mon=mon+1
  endwhile
  return, [mon,dum]
end
