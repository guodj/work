function mkupper, string

  temp = byte(string)

  loc = where((temp ge 97) and (temp le 122), count)

  if count ne 0 then temp(loc) = temp(loc)-32

  return, string(temp)

end

