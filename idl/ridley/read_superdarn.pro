
pro read_superdarn, file, time, potential, lats, mlts

  iStartTime = intarr(6)
  iEndTime = intarr(6)

  nPtsSingle = 0L
  nPtsTotal = 0L

  nPtsMax = 1440L*5

  dLat = 1.0
  dLon = 2.0

  MinLat = 60.0

  nLats = (90.0-MinLat-dLat)/dLat + 1
  nLons = 360.0/dLon

  nPtsExpected = nLons*nLats

  potential = fltarr(nLons,nLats)
  lats      = fltarr(nLons,nLats)
  mlts      = fltarr(nLons,nLats)

  potential_rot = fltarr(nLons+1,nLats+1)
  lats_rot      = fltarr(nLons+1,nLats+1)
  mlts_rot      = fltarr(nLons+1,nLats+1)

  potential_out = fltarr(nPtsMax,nLons+1,nLats+1)
  lats_out = fltarr(nLons+1,nLats+1)
  mlts_out = fltarr(nLons+1,nLats+1)

  time = dblarr(nPtsMax)

  for i=0,nLats do mlts_rot(*,i) = findgen(nLons+1)/nLons*24.0

  close,1
  openr,1,file

  tmp = fltarr(7)

  while not eof(1) do begin
     
     readf,1,iStartTime,iEndTime
     readf,1,nPtsSingle

;     print, iStartTime
;     print, iEndTime
;     print, nPtsSingle

     c_a_to_r, iStartTime, sTime
     c_a_to_r, iEndTime, eTime
     print, nPtsTotal
     time(nPtsTotal) = (stime+etime)/2.0

;     if (nPtsSingle ne nPtsExpected) then begin
;        print, 'nPtsSingle ne nPtsExpected. Please see the code!'
;        print, 'Single/Expected : ',nPtsSingle, nPtsExpected
;     endif

     potential(*,*) = 0.0

     for iPts = 0,nPtsSingle-1 do begin

        readf,1,tmp
        iLat = (tmp(4)-MinLat)/dLat
        iLon = tmp(5)/dLon

        potential(iLon,iLat) = tmp(0)
        lats(iLon,iLat) = tmp(4)
        mlts(iLon,iLat) = tmp(6)

     endfor

     m = mlts(*,0)
     l = where(m eq min(m))
     l = l(0)
     dmlt = m((l+1) mod nLons)-m(l)

     ii = fix(l-1+nLons) mod nLons
     jj = (ii + 1) mod nLons
     dx = (24.0-mlts(ii)) / dmlt

     mltnew = (1.0-dx)*(mlts(ii,0)-24.0) + dx*mlts(jj,0)
     
;     print, ii, jj, dx, mlts(ii), mltnew
;     stop

     potential_rot(0,0:nLats-1) = (1.0-dx)*potential(ii,*) + dx*potential(jj,*)
     lats_rot(0,0:nLats-1) = (1.0-dx)*lats(ii,*) + dx*lats(jj,*)
        
     for i=1,nLons-1 do begin

        ii = fix(mlts_rot(i,0)/dmlt + l - 1) mod nLons
        jj = (ii + 1) mod nLons
        dx = (mlts_rot(i,0) - mlts(ii)) / dmlt

        mltnew = (1.0-dx)*mlts(ii,0) + dx*mlts(jj,0)

        potential_rot(i,0:nLats-1) = (1.0-dx)*potential(ii,*) + dx*potential(jj,*)
        lats_rot(i,0:nLats-1) = (1.0-dx)*lats(ii,*) + dx*lats(jj,*)
        
     endfor

     potential_rot(nLons,*) = potential_rot(0,*)
     lats_rot(nLons,*) = lats_rot(0,*)
     mlts_rot(nLons,*) = mlts_rot(0,*)

     potential_rot(1,*) = (potential_rot(0,*) + potential_rot(2,*))/2.0
     potential_rot(*,nLats) = mean(potential_rot(*,nLats-1))
     lats_rot(*,nLats) = 90.0

     ;x = (90.0-lats_rot)*cos(mlts_rot*!pi/12.0-!pi/2) 
     ;y = (90.0-lats_rot)*sin(mlts_rot*!pi/12.0-!pi/2)
     ;contour, potential_rot/1000.0, x, y, /follow, nlevels = 10

     potential_out(nPtsTotal,*,*) = potential_rot

     nPtsTotal++

  endwhile

  close,1

  time = time(0:nPtsTotal-1)
  potential = potential_out(0:nPtsTotal-1,*,*)

  lats = lats_rot
  mlts = mlts_rot

end
