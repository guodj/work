function sgn, num
  val = num
  loc = where(num ne 0.0,count)
  if count gt 0 then val(loc) = num(loc)/abs(num(loc))
  loc = where(num eq 0.0,count)
  if count gt 0 then val(loc) = 1.0
  return, val
end
