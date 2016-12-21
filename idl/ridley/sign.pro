
function sign, array

  s = array*0.0 + 1.0
  loc = where(array lt 0.0,count)
  if (count gt 0) then s(loc) = -1.0
  return, s

end 
