
pro iono_read_file_tec, filename, nvars, nlats, nlons, vars, time, data

  line = ""
  nVars = 0
  nLats = 0
  nLons = 0

  close,1
  openr,1,filename

  done = 0

  readf,1,line

  l = strpos(line,'-')
  year   = fix(strmid(line,l-4,4))
  month  = fix(strmid(line,l+1,2))
  day    = fix(strmid(line,l+4,2))
  hour   = fix(strmid(line,l+7,2))
  minute = fix(strmid(line,l+10,2))
  second = fix(strmid(line,l+13,2))

  print, year, month, day, hour, minute, second

  itime = [year, month, day, hour, minute, second]
  c_a_to_r, itime, time

  readf,1,line

  print,line
  l = strpos(line,'=')
  line = strmid(line,l+1,strlen(line)-l)

  vars = strsplit(line,',',/extract)

  while (strpos(line,'ZONE') lt 0) do begin

     readf,1,line
     if (strpos(line,'ZONE') lt 0) then begin
        subvar = strsplit(line,',',/extract)
        vars = [vars,subvar]
     endif

  endwhile

  nVars = n_elements(vars)
  for iVar = 0,nVars-1 do begin
     v = vars(iVar)
     l = strpos(v,'"')
     v = strmid(v,l+1,strlen(v)-l)
     l = strpos(v,'"')
     v = strmid(v,0,l)
     vars(iVar) = v
  endfor

  readf,1,line

  l = strpos(line,'I=')
  nLats = fix(strmid(line,l+2,strlen(line)))
  l = strpos(line,'J=')
  nLons = fix(strmid(line,l+2,strlen(line)))

  data = fltarr(2,nvars,nlons,nlats)

  tmp = fltarr(nVars)

  for iLon=0,nLons-1 do for iLat=0,nLats-1 do begin
     readf,1,tmp
     data(0,*,iLon,iLat) = tmp
  endfor

  readf,1,line
  readf,1,line

  for iLon=0,nLons-1 do for iLat=0,nLats-1 do begin
     readf,1,tmp
     data(1,*,iLon,iLat) = tmp
  endfor

  close,1

end

