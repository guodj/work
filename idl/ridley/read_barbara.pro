pro read_barbara, filein, ltpos, lnpos, data, time, date, lats, lons,	$
    type, spd, by, bz

  close,1
  openr,1,filein
 
  line = ''
  readf,1, line
  start = strpos(line, "AMIE")-20
  lons = fix(strmid(line,start+16,3))+1
  lats = fix(strmid(line,start+12,3))+1
  count = fix(strmid(line,start,4))
 
  format = '('+tostr(lons)+'E12.4,A5)'

  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  lts = fltarr(lats)
  fakedata = fltarr(1,lons)
  date = strarr(count)
  time = strarr(count)
  data = fltarr(count,lons,lats)
  by   = strarr(count)
  bz   = strarr(count)
  spd  = strarr(count)
  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'
 
  done = 0
  while not done do begin
 
    readf,1, line
    if strpos(line,"ntime iyear") gt 0 then done = 1
 
  endwhile
 
  type = ''

  for i=0,count-1 do begin
 
    year = strmid(line,5,4)
    month = fix(strmid(line,10,2))
    date(i) = strmid(mon,(month-1)*3,3)+' '+strmid(line,13,2)+', '+year
    h = fix(strmid(line,16,2))
    m = fix(strmid(line,19,2))
    time(i) = chopr('0'+tostr(h),2)+':'+chopr('0'+tostr(m),2)+' UT'
 
    readf,1, line
    spd(i) = 'Speed : '+strcompress(strmid(line, 0,6),/remove)
    by(i)  = 'By : '+strcompress(strmid(line,12,6),/remove)+' nT'
    bz(i)  = 'Bz : '+strcompress(strmid(line,18,6),/remove)+' nT'
    readf,1, line

    for j=0,lats-1 do begin
 
      readf,1, line
      lts(j) = float(strmid(line,5,5))
      readf,1, format=format,fakedata, type
      data(i,0:lons-1,j) = fakedata(0,0:lons-1)
 
    endfor
 
    if i lt count-1 then readf,1, line
 
  endfor

  close,1

  for i=0,lons-1 do for j=0,lats-1 do begin
    ltpos(i,j) = 90.0 - lts(j)
    lnpos(i,j) = float(i)*360.0/float(lons-1)
  endfor
 
  loc = where(lnpos lt 0.0, count)
 
  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  type = strmid(type,1,strlen(type)-1)

  if (strpos(type,'pot') gt -1) or (strpos(type,'efpo') gt -1) then begin
    data = data/1000.0
    type = 'Electric Potential'
  endif
  if strpos(type,'sigp') gt -1 then type = 'Pedersen Conductance'
  if strpos(type,'sigh') gt -1 then type = 'Hall Conductance'
  if strpos(type,'dfac') gt -1 then type = 'Field Aligned Current'
  if strpos(type,'sjht') gt -1 then type = 'Simple Joule Heating'
  if strpos(type,'efld') gt -1 then begin
    type = '(Eastward) E-Field (mV/m)'
    data = data*1000
  endif
  if strpos(type,'nfld') gt -1 then begin
    type = '(Northward) E-Field (mV/m)'
    data = data*1000
  endif
  return

end
