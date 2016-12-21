
pro read_champ, massfile, posfile, time, position, MassDensity, dir = dir, $
                iYear = iYear, iDay = iDay

  if (n_elements(dir) eq 0) then dir = '.'

  if (n_elements(iYear) eq 0) then iYear = fix(strmid(posfile,0,4))
  if (n_elements(iDay) eq 0) then iDay  = fix(strmid(posfile,5,3))

  iTime = [iYear, 1, iDay, 0,0,0]

  c_a_to_r, iTime, BaseTime

  close,1
  openr,1, dir+'/'+posfile

  line = ''
  readf,1,line

  hour = 0.0
  lat  = 0.0
  lon  = 0.0
  alt  = 0.0
  mass = 0.0

  nPts = 0
  position = fltarr(3,30000)
  time     = dblarr(30000)
  day      = 0.0
  oldhour  = -1.0

  while not eof(1) do begin

      readf,1, hour, lat, lon, alt
      lon = (lon + 360.0) mod 360.0

      if (hour lt oldhour) then day = day + 1.0
      oldhour = hour

      rTime = BaseTime + day*24.0*3600.0 + hour*3600.0
      c_r_to_a, iTime, rTime
      position(0,nPts) = lon
      position(1,nPts) = lat
      position(2,nPts) = alt
      time(nPts)       = rTime
      nPts             = nPts + 1

  endwhile

  close,1

  nPts = nPts - 1

  MassDensity = fltarr(nPts+1)

  openr,1,dir+'/'+MassFile

  for iPt=0,nPts-1 do begin

      readf,1,hour,lat,lon,mass
      MassDensity(iPt) = mass*1.0e-12

  endfor

  close,1

  MassDensity = MassDensity(0:nPts-1)
  Time        = Time(0:nPts-1)
  Position    = Position(*, 0:nPts-1)

end

