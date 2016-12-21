
pro read_goce, dir, data, starttime, endtime

  c_r_to_a, itime, starttime

  year  = tostr(itime(0),4)
  month = tostr(itime(1),2)

  file = dir+'/goce_denswind_v1_4_'+year+'-'+month+'.txt'

  spawn, 'wc '+file, wc
  nPtsMax = long(wc(0))

  time = dblarr(nPtsMax)
  alt  = dblarr(nPtsMax)
  lat  = fltarr(nPtsMax)
  lon  = fltarr(nPtsMax)
  density = fltarr(nPtsMax)
  east = fltarr(nPtsMax)
  north = fltarr(nPtsMax)
  vertical = fltarr(nPtsMax)
  denerror = fltarr(nPtsMax)
  winderror = fltarr(nPtsMax)


  close,1
  openr,1,file

  line = ''

  nPts = 0L

  IsDone = 0

  while not IsDone do begin

     readf,1,line

     if eof(1) then IsDone = 1 $
     else begin

        if (strpos(line,'#') eq -1) then begin

           members = strsplit(line,' ',/extract)
           ymd = fix(strsplit(members(0),'-',/extract))
           hms = fix(strsplit(members(1),':',/extract))
           itime = [ymd,hms]
           c_a_to_r, itime, rtime
           if (rtime ge starttime and rtime lt endtime) then begin

              alt(nPts) = double(members(3))
              lat(nPts) = float(members(5))
              lon(nPts) = float(members(4))
              density(nPts) = float(members(8))
              east(nPts) = float(members(9))
              north(nPts) = float(members(10))
              vertical(nPts)  = float(members(11))
              denerror(nPts)  = float(members(12))
              winderror(nPts) = float(members(13))
              time(nPts) = rtime
              
              nPts++

           endif

           if (rTime gt endtime) then IsDone = 1

        endif

     endelse

  endwhile

  close,1

  if (nPts eq 0) then nPts = 1

  alt = alt(0:nPts-1)
  lat = lat(0:nPts-1)
  lon = (lon(0:nPts-1)+360.0) mod 360.0
  density = density(0:nPts-1)
  east = east(0:nPts-1)
  north = north(0:nPts-1)
  vertical = vertical(0:nPts-1)
  denerror = denerror(0:nPts-1)
  winderror = winderror(0:nPts-1)

  time = time(0:nPts-1)

  data = {nPts: nPts, $
          alt: alt, $
          lat: lat, $
          lon: lon, $
          density: density, $
          east: east, $
          north: north, $
          vertical: vertical, $
          denerror: denerror, $
          winderror: winderror, $
          time: time}

end
