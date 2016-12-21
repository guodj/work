
pro thermo_sza, lats, lons, CurrentTime, LocalTime, sza

  itime = [1999, 3, 21, 0, 0, 0]
  c_a_to_r, itime, VernalTime

  SecondsPerYear = 365.25*24.0*60.0*60.0
  Tilt = 23.5*!dtor

  OrbitAngle = 2*!pi*(CurrentTime-VernalTime)/SecondsPerYear
  SunDeclination = atan(tan(Tilt)*sin(OrbitAngle))

  SinDec = sin(SunDeclination)
  CosDec = cos(SunDeclination)

  Ut = CurrentTime mod (24.0*60.0*60.0)

  LocalTime = (Ut/3600.0 + lons * 24.0 / 360.0) mod 24.0

  sza = acos(SinDec*sin(lats*!dtor) + $
             CosDec*cos(lats*!dtor)*cos(!pi*(LocalTime-12.0)/12.0))/!dtor

end
