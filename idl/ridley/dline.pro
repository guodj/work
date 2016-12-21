pro dline,line

  for i=0,strlen(line)-1 do begin
    print,i,' - ',strmid(line,i,1)
  endfor

  return

end
