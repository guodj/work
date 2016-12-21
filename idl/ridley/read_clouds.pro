
pro read_clouds, file, time, data, location, station, statlist

  line = ''

  nPtsMax = 31*25
  nPts = 0

  time     = dblarr(nPtsMax)
  data     = fltarr(nPtsMax)
  station  = fltarr(nPtsMax)
  statlist = strarr(5)

  openr,1,file
  
  IsDone = 0
  itime = intarr(6)
  while (not IsDone) do begin

     readf,1,line
     if (strpos(line,'Latitude') gt -1) then readf,1,lat
     if (strpos(line,'Longitude') gt -1) then readf,1,lon
     if (strpos(line,'Stations Used') gt -1) then begin
        for iStat = 0, 4 do begin
           readf,1,line
           statlist(iStat) = line
        endfor
     endif
     if (strpos(line,'CLOUDCOVER') gt -1) then begin
        IsDone = 1
        tmp = intarr(6)
        while not eof(1) do begin
           readf,1,tmp
           itime(0:3) = tmp(0:3)
           c_a_to_r, itime, rtime
           DoAddPoint = 1
           if (nPts gt 0) then begin
              if (rtime eq time(nPts-1)) then begin
                 data(nPts-1) = (data(nPts-1)+tmp(5))/2
                 DoAddPoint = 0
              endif
           endif
           if (DoAddPoint) then begin
              time(nPts) = rtime
              station(nPts) = tmp(4)
              data(nPts) = tmp(5)
              nPts++
           endif
        endwhile
     endif

  endwhile

  close,1

  time = time(0:nPts-1)
  data = data(0:nPts-1)
  station = station(0:nPts-1)

end
