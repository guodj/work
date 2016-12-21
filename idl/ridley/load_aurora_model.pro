
pro read_value, unit, nLats, nMlts, nLevels, value

  line = ''
  readf,1,line
  tmp = fltarr((nMlts-1)/2)
  for iLevel=0,nLevels-1 do begin
     for iLat = 0,nLats-1 do begin
        readf,1,tmp
        Value(iLevel,0:(nMlts-1)/2-1,iLat) = tmp
        readf,1,tmp
        Value(iLevel,(nMlts-1)/2:nMlts-2,iLat) = tmp
        Value(iLevel,nMlts-1,iLat) = Value(iLevel,0,iLat)
     endfor
  endfor
  Value = Value * 0.001

end

pro load_aurora_model, aurora


  file = '/raid3/Data/HemisphericPower/hpke.noaa'

  openr,1,file

  nLats = 21
  nMlts = 31
  nLevels = 10

  line = ''
  for i=1,4 do readf,1,line

  Hall  = fltarr(nLevels, nMlts, nLats)
  Ped   = fltarr(nLevels, nMlts, nLats)
  AveE  = fltarr(nLevels, nMlts, nLats)
  EFlux = fltarr(nLevels, nMlts, nLats)

  read_value, unit, nLats, nMlts, nLevels, Hall
  read_value, unit, nLats, nMlts, nLevels, Ped
  read_value, unit, nLats, nMlts, nLevels, AveE
  read_value, unit, nLats, nMlts, nLevels, EFlux

  close,1

  Lats = fltarr(nMlts, nLats)
  Mlts = fltarr(nMlts, nLats)

  for iLat = 0,nLats-1 do Mlts(*,iLat) = findgen(nMlts)/(nMlts-1)*24.0
  for iMlt = 0,nMlts-1 do Lats(iMlt,*) = findgen(nLats)/(nLats-1)*40.0+50.0

  x = (90.0-Lats)*cos(Mlts*!pi/12-!pi/2)
  y = (90.0-Lats)*sin(Mlts*!pi/12-!pi/2)

  area = 120000.0^2 * (360.0/(nMlts-1)) * (40.0/(nLats-1)) * cos(Lats*!dtor)

  hp = fltarr(nLevels)
  for iLevel = 0,nLevels-1 do begin
     hp(iLevel) = total(area*eflux(iLevel,*,*)/1000.0)/1e9
  endfor

  aurora = { $
           nLats : nLats, $
           nMlts : nMlts, $
           nLevels : nLevels, $
           Hall : Hall, $
           Ped  : Ped, $
           EFlux : EFlux, $
           AveE  : AveE, $
           Mlts  : Mlts, $
           Lats  : Lats, $
           Area  : Area, $
           x:x, y:y, $
           hp : hp}

end



