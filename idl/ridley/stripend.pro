function stripend, line

  strl = strlen(line)

  done = 0

  while not done do begin
    if strmid(line,strl-1,1) eq ' ' then strl = strl - 1 else done = 1
  endwhile

  outl = strmid(line,0,strl)

  return, outl

end