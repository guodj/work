function get_month, month

  m = 'JanFebMarAprMayJunJulAugSepOctNovDec'

  r = strmid(m,(month-1)*3,3)

  return, r

end
