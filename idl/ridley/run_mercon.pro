PRO run_mercon

  print,'WARNING: You must be in the directory with the files created by'
  print,'remmean.f. If you are not, break the program, exit IDL, and move'
  print,'the correct directory.'

  print, ' '
  print, 'Enter center of file name (ex 9105 for sss9105.remca)'
  print, 'Enter return for *.remca) :'
  middle = ''
  read, middle
  if strlen(middle) gt 0 then fname = '*'+middle+'.remca' $
  else fname = '*.remca'

  print, 'Enter filename with mlats in it :'
  latfile = ''
  read, latfile

  close,1
  openr,1,latfile
  line = ''
  readf,1,line
  ntpd = 1440/fix(strmid(line,14,3))

  print, 'Enter dst file (Downloaded html file from Kyoto):'
  dstfile = ''
  read,dstfile

  print, 'Enter quiet day and quiet day hour to start at (day,hour) :'
  day = 0
  hour = 0
  read, day, hour
  if hour eq 0 then hour = 1

  print, 'Enter start day, end day :'
  sday = 0
  eday = 0
  read, sday, eday

  nday = eday - sday + 1

  qdst = intarr(1,24)
  adst = intarr(nday,24)

  list = findfile(fname)

  nstat = n_elements(list)

  close,1

  col = indgen(24)*4 + 3
  col(8:23) = col(8:23) + 1
  col(16:23) = col(16:23) + 1

  openr,1,dstfile
  line = ''
  done = 0
  while (not done) do begin
    readf,1,line
    if strpos(line,"DAY") eq 0 then done = 1
  endwhile
  done = 0
  while (not done) do begin
    readf,1,line
    if fix(strmid(line,0,2)) eq day then done = 1
  endwhile
  for i=hour,24 do qdst(i-1) = fix(strmid(line,col(i-1),4))
  if (hour gt 1) then begin
    readf,1,line
    if (strlen(line) lt 40) then readf, 1, line 
    for i=1,hour-1 do qdst(i-1) = fix(strmid(line,col(i-1),4))
  endif
  close,1

  openr,1,dstfile
  line = ''
  done = 0
  while (not done) do begin
    readf,1,line
    if strpos(line,"DAY") eq 0 then done = 1
  endwhile
  done = 0
  while (not done) do begin
    readf,1,line
    if fix(strmid(line,0,2)) eq sday then done = 1
  endwhile
  for j=sday,eday do begin
    for i=1,24 do adst(j-sday,i-1) = fix(strmid(line,col(i-1),4))
    readf,1,line
    if (strlen(line) lt 40) then readf, 1, line 
  endfor
  close,1

  if (strlen(middle) gt 0) then tfile = 'dst'+middle+'.amiein' $
  else tfile = 'dst.amiein'
  openw,1,tfile

  for i=sday,eday do begin

    ds = strcompress(string(i))

    for j=0,5 do begin

      h1 = ' '+strcompress(string(j*4),/remove_all)
      h2 = ' '+strcompress(string(j*4+1),/remove_all)
      h3 = ' '+strcompress(string(j*4+2),/remove_all)
      h4 = ' '+strcompress(string(j*4+3),/remove_all)
      h1 = strmid(h1, strlen(h1)-2,2)
      h2 = strmid(h2, strlen(h2)-2,2)
      h3 = strmid(h3, strlen(h3)-2,2)
      h4 = strmid(h4, strlen(h4)-2,2)

      d1 = '     '+strcompress(string(adst(i-sday,j*4)),/remove_all)
      d2 = '     '+strcompress(string(adst(i-sday,j*4+1)),/remove_all)
      d3 = '     '+strcompress(string(adst(i-sday,j*4+2)),/remove_all)
      d4 = '     '+strcompress(string(adst(i-sday,j*4+3)),/remove_all)
      d1 = strmid(d1, strlen(d1)-3,3)+'.0'
      d2 = strmid(d2, strlen(d2)-3,3)+'.0'
      d3 = strmid(d3, strlen(d3)-3,3)+'.0'
      d4 = strmid(d4, strlen(d4)-3,3)+'.0'
      if adst(i-sday,j*4) le -100.0 then 				$
        d1 = strcompress(string(fix(adst(i-sday,j*4)))+'.',/remove_all)
      if adst(i-sday,j*4+1) le -100.0 then 				$
        d2 = strcompress(string(fix(adst(i-sday,j*4+1)))+'.',/remove_all)
      if adst(i-sday,j*4+2) le -100.0 then 				$
        d3 = strcompress(string(fix(adst(i-sday,j*4+2)))+'.',/remove_all)
      if adst(i-sday,j*4+3) le -100.0 then 				$
        d4 = strcompress(string(fix(adst(i-sday,j*4+3)))+'.',/remove_all)

      printf,1,' '+ds+','+h1+',30,'+d1+','+ds+','+h2+',30,'+d2+','+	$
	       ds+','+h3+',30,'+d3+','+ds+','+h4+',30,'+d4+','

    endfor

  endfor

  close,1

  openw,1,'temp.input'

  printf,1,strcompress(string(nstat))
  printf,1,strcompress(string(nday))
  printf,1,strcompress(string(sday))
  printf,1,strcompress(string(ntpd))

  for i=0,nstat-1 do printf,1,list(i)

  for i=0,23 do printf,1, strcompress(string(qdst(0,i)))

  for j=0,nday-1 do for i=0,23 do printf,1, strcompress(string(adst(j,i)))

  printf,1,latfile
  if (strlen(middle) gt 0) then $
    tfile = 'mag'+middle+strcompress(string(sday),/remove_all)+'.final' $
  else tfile = 'mag'+strcompress(string(sday),/remove_all)+'.final'
  printf,1,tfile

  if (strlen(middle) gt 0) then $
    tfile = 'ae'+middle+strcompress(string(sday),/remove_all)+'.final' $
  else tfile = 'ae'+strcompress(string(sday),/remove_all)+'.final'
  printf,1,tfile

  close,1

  print,'To run mercon now, simply type "mercon < temp.input" at the prompt'

  return

end
