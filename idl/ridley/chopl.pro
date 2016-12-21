function chopl, svalue, n
  if strlen(svalue) lt n then n = strlen(svalue)
  return, strmid(svalue, 0, n)
end
