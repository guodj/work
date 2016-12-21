;-------------------------------------------------------------------------
;  This procedure finds the index of the magnetometer station name
;  if it is already in the list, else it adds it to the list and
;  returns the last index.
;-------------------------------------------------------------------------

pro findmag, MagList, nMagsMax, stat, iMag

  iMag = 0
  stat = strcompress(stat,/remove)
  while (iMag lt nMagsMax and $
            strlen(MagList(iMag)) gt 0 and $
            strpos(MagList(iMag), stat) lt 0) do iMag = iMag + 1
  if (strlen(MagList(iMag)) eq 0) then MagList(iMag) = stat

end


pro read_amie_data, File, RawMagData, AMIEMagData, MagFlag, MinFlag, $
                    MagTime, MagList, nMags, nMagsTotal

  nMagsMax    = 200
  nTimesMax   = 1440*10
  RawMagData  = fltarr(nMagsMax,3,nTimesMax) - 9999.0
  AMIEMagData = fltarr(nMagsMax,3,nTimesMax) - 9999.0
  MagFlag     = intarr(nMagsMax,3,nTimesMax)
  MinFlag     = intarr(nTimesMax)
  MagTime     = dblarr(nTimesMax)
  MagList     = strarr(nMagsMax)
  nMags       = intarr(nTimesMax)
  
  close,1
  line = ""
  iTimeArray = intarr(7)

  mlat   = 0.0
  mlt    = 0.0
  h      = 0.0
  amieh  = 0.0
  d      = 0.0
  amied  = 0.0
  z      = 0.0
  amiez  = 0.0
  iFlagH = 0
  iFlagD = 0
  iFlagZ = 0
  stat   = ""

  fl = findfile(File)
  if (strlen(fl(0)) eq 0) then return

  ; ------------------------------------------------------------------------
  ; This section of code reads in the magnetometer data from the _data files
  ; ------------------------------------------------------------------------

  ; we need to figure out whether this is an old-style summary file, or
  ; a new file.

  spawn, "grep MAGFLAG "+File, list
  if (strlen(list(0)) eq 0) then IsOldFile = 1 else IsOldFile = 0

  openr,1, File

  iTime = -1

  IsDone = 0

  while not IsDone do begin

      readf,1,line

      ; Read in the time

      if (strpos(line,"#TIME") gt -1) then begin

          ; We up the iteration number (iTime) when the time variable
          ; is read in, since this is the first thing output by AMIE
          iTime = iTime + 1
          if (iTime ge nTimesMax) then begin
              print, "iTime > nTimesMax : ", iTime, nTimesMax
              stop
          endif

          readf,1,iTimeArray
          ; Convert the time to a real number, and store it
          c_a_to_r, iTimeArray(0:5), rTime
          MagTime(iTime) = rTime
          line = ""

      endif

      ; Read in the magnetometer data

      if (strpos(line,"#MAGNETOMETERS") gt -1) then begin

          readf,1,n
          nMags(iTime) = n

          for i = 0, n-1 do begin
              if (IsOldFile) then begin
                  readf,1,mlat,mlt,h,amieh,d,amied,z,amiez
                  iFlagH = 0
                  iFlagD = 0
                  iFlagZ = 0
                  stat = tostr(i)
              endif else begin
                  readf,1,mlat,mlt,h,amieh,d,amied,z,amiez, $
                    iFlagH,iFlagD,iFlagZ,stat
              endelse
              findmag, MagList, nMagsMax, stat, iMag
              RawMagData(iMag,0,iTime) = h
              RawMagData(iMag,1,iTime) = d
              RawMagData(iMag,2,iTime) = z
              AMIEMagData(iMag,0,iTime) = amieh
              AMIEMagData(iMag,1,iTime) = amied
              AMIEMagData(iMag,2,iTime) = amiez
              MagFlag(iMag,0,iTime) = iFlagH
              MagFlag(iMag,1,iTime) = iFlagD
              MagFlag(iMag,2,iTime) = iFlagZ
          endfor

          if (IsOldFile) then begin
              nM = 0
              for i=0,nMags(iTime)-1 do begin
                  if (max(abs(RawMagData(i,*,iTime))) lt 9998.0) then nM = nM+1
              endfor
              nMags(iTime) = nM
          endif

          line = ""
      endif

      ; Read in the Flags - these are only there if they are 1 or more.

      if (strpos(line,"#MAGFLAGS") gt -1) then begin
          readf,1, n
          readf,1, iMinErrorFlag
          MinFlag(iTime) = iMinErrorFlag
          ;for i = 0, n-1 do begin
          ;    readf,1,format="(a3,3i3)",stat,iFlagH,iFlagD,iFlagZ
          ;    findmag, MagList, nMagsMax, stat, iMag
          ;    MagFlag(iMag,0,iFile,iTime) = iFlagH
          ;    MagFlag(iMag,1,iFile,iTime) = iFlagD
          ;    MagFlag(iMag,2,iFile,iTime) = iFlagZ
          ;endfor

          nM = 0
          for i=0,nMags(iTime)-1 do begin
              if (max(MagFlag(i,*,iTime)) le MinFlag(iTime) and $
                  max(abs(RawMagData(i,*,iTime))) lt 9998.0) then nM = nM+1
          endfor
          nMags(iTime) = nM

          line = ""
      endif

      if (eof(1)) then IsDone = 1

  endwhile

  close,1

  nMagsTotal = 0
  while (strlen(MagList(nMagsTotal)) gt 0) do nMagsTotal = nMagsTotal+1


end
