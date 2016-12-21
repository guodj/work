function pad, stringval, numtopad

  ds = stringval
  for i=1,numtopad do ds = ds + ' '
  ds = strmid(ds,0,numtopad)

  return, ds

end
