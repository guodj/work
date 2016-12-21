
pro read_dmsp_webfile, dmspfile, satellite, nPts, $
                       time, Quality, GeoPos, MagPos, Velocity, $
                       Sigma, NDensity, Frac, Temperature, Points

  line        = ''
  tmp         = dblarr(22)
  nPtsMax     = 10000
  time        = dblarr(nPtsMax)
  quality     = intarr(2,nPtsMax)
  GeoPos      = fltarr(3,nPtsMax)
  MagPos      = fltarr(2,nPtsMax)
  Velocity    = fltarr(3,nPtsMax)
  Sigma       = fltarr(3,nPtsMax)
  NDensity    = fltarr(nPtsMax)
  Frac        = fltarr(3,nPtsMax)
  Temperature = fltarr(2,nPtsMax)
  Points      = intarr(2,nPtsMax)
  iTime       = intarr(6)

  close, 1
  openr, 1, dmspfile

  readf,1,line
  satellite = strmid(line,0,3)
  readf,1,line
  readf,1,line

  nPts = 0
  while not eof(1) do begin

      readf,1,tmp

      iTime(0) = (tmp(0)/1000.0) mod 100
      iTime(1) = 1
      iTime(2) = tmp(0) mod 1000

      iHour = fix(tmp(1)/3600.0)
      Left  = tmp(1) - iHour*3600.0

      iTime(3) = iHour
      iTime(4) = 0
      iTime(5) = Left

      c_a_to_r, iTime, rTime
      time(nPts) = rTime

      Quality(*,nPts)     = tmp(2:3)

      GeoPos(0,nPts)      = tmp(6)
      GeoPos(1,nPts)      = tmp(5)
      GeoPos(2,nPts)      = tmp(4)

      MagPos(0,nPts)      = tmp(8)
      MagPos(1,nPts)      = tmp(7)

      Velocity(*,nPts)    = tmp(9:11)
      Sigma(*,nPts)       = tmp(12:14)

      NDensity(nPts)      = tmp(15)
      Frac(*,nPts)        = tmp(16:18)

      Temperature(*,nPts) = tmp(19:20)

      Points(nPts)        = tmp(21)

      nPts = nPts + 1

  endwhile

  close,1

  nPts = nPts - 1

  time        = time(0:nPts-1)
  quality     = quality(*,0:nPts-1)
  GeoPos      = GeoPos(*,0:nPts-1)
  MagPos      = MagPos(*,0:nPts-1)
  Velocity    = Velocity(*,0:nPts-1)
  Sigma       = Sigma(*,0:nPts-1)
  NDensity    = NDensity(0:nPts-1)
  Frac        = Frac(*,0:nPts-1)
  Temperature = Temperature(*,0:nPts-1)
  Points      = Points(0:nPts-1)

end
