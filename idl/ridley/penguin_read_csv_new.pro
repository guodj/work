
pro penguin_read_csv_new, filename, Vars, time, data, error

  spawn, 'wc '+filename, wc
  nLines = long(wc(0))-1

  close,1
  openr,1,filename

  line = ''
  readf,1,line
  vars = strsplit(line,',',/extract)

  nVars = n_elements(vars)
  if (nVars lt 3 or nVars gt 4) then begin
     print, 'In penguin_read_csv_new:'
     print, 'nVars ne 3! line = ',line
     stop
  endif else nVars = 3

  data = fltarr(nLines,3)
  time = dblarr(nLines)
  error = intarr(nLines)

  t = strsplit(filename,'_',/extract)
  nT = n_elements(t)
  iTime = fix(t(nT-6:nT-1))
  c_a_to_r, iTime, basetime

  for i=0L, nLines-1 do begin
     readf,1,line
     d = strsplit(line,',',/extract)
     for j=0,2 do data(i,j) = float(d(j))
     error(i) = fix(d(3))
     time(i) = basetime + double(i)
  endfor

  close,1

end
