

pro fpi_read_makela, file, fpi, alldata

  spawn,'wc '+file, wc
  nLines = fix(wc(0))-3

  close,1
  openr,1,file
  line = ''
  readf,1,line
  readf,1,line

  readf,1,line
  vars = strsplit(line,' ',/extract)

  ;display, vars

  nVars = n_elements(vars)-4
  vars = vars(1:nVars-4)

  fpidata = fltarr(nVars,nLines)
  fpitime = dblarr(nLines)

  for iLine = 0, nLines-1 do begin

     readf,1,line
     d = strsplit(line,' ',/extract)
     date = fix(strsplit(d(0),'-',/extract))
     time = fix(strsplit(d(1),':',/extract))
     itime = [date,time]
     c_a_to_r, itime, rtime
     fpitime(iLine) = rtime

     fpidata(*,iLine) = float(d(2:nVars+1))

  endfor

  close,1

  fpiAz = round(reform(fpidata(0,*))/90.0)
  fpiZe = round(reform(fpidata(1,*))/45.0)
  fpiT  = reform(fpidata(2,*))
;  fpiTs = reform(fpidata(3,*))
  fpiV  = reform(fpidata(4,*))
;  fpiVs = reform(fpidata(5,*))
  fpiTflag = reform(fpidata(15,*))
  fpiVflag = reform(fpidata(15,*))

  alldata = { $
            time: fpitime, $
            azimuth:reform(fpidata(0,*)), $
            zenith: reform(fpidata(1,*)), $
            temperature: fpiT, $
            velocity: fpiV, $
            FlagT: fpiTflag, $
            FlagV: fpiVflag}

  l = where(fpiZe eq 0,c)
  if (c gt 0) then begin
     zenith = fpiV(l)
     zenithTime = fpitime(l)
     zenithTemp = fpiT(l)
     zenithTflag = fpiTflag(l)
     zenithVflag = fpiVflag(l)
  endif else begin
     zenith = [0.0]
     zenithTime = [min(fpitime(0))]
     zenithTemp = [0.0]
     zenithTflag = [0.0]
     zenithVflag = [0.0]
  endelse
  
; East may be az = 1; ze = 1
; West may be az = 3; ze = 1
; North may be az = 0; ze = 1
; South may be az = 2; ze = 1

l = where(abs(fpiAz) eq 0 and fpiZe eq 1,c)
if (c gt 0) then begin
   los = fpiV(l)
   NorthTime = fpitime(l)
   w = interpolate_mine(NorthTime, zenith, zenithTime)
   ze = !pi/4
   north = (los - w*cos(ze))/sin(ze)
   northTemp  = fpiT(l)
   northTflag = fpiTflag(l)
   northVflag = fpiVflag(l)
endif else begin
   north = [0.0]
   northTime = [min(fpitime(0))]
   northTemp = [0.0]
   northTflag = [0.0]
   northVflag = [0.0]
endelse

l = where(abs(fpiAz) eq 2 and fpiZe eq 1,c)
if (c gt 0) then begin
   los = fpiV(l)
   SouthTime = fpitime(l)
   w = interpolate_mine(SouthTime, zenith, zenithTime)
   ze = !pi/4
   south = - (los - w*cos(ze))/sin(ze)
   southTemp  = fpiT(l)
   southTflag = fpiTflag(l)
   southVflag = fpiVflag(l)
endif else begin
   south = [0.0]
   southTime = [min(fpitime(0))]
   southTemp = [0.0]
   southTflag = [0.0]
   southVflag = [0.0]
endelse

l = where(abs(fpiAz) eq 1 and fpiZe eq 1,c)
if (c gt 0) then begin
   los = fpiV(l)
   EastTime = fpitime(l)
   w = interpolate_mine(EastTime, zenith, zenithTime)
   ze = !pi/4
   east = (los - w*cos(ze))/sin(ze)
   eastTemp  = fpiT(l)
   eastTflag = fpiTflag(l)
   eastVflag = fpiVflag(l)
endif else begin
   east = [0.0]
   eastTime = [min(fpitime(0))]
   eastTemp = [0.0]
   eastTflag = [0.0]
   eastVflag = [0.0]
endelse

l = where(abs(fpiAz) eq 3 and fpiZe eq 1,c)
if (c gt 0) then begin
   los = fpiV(l)
   WestTime = fpitime(l)
   w = interpolate_mine(WestTime, zenith, zenithTime)
   ze = !pi/4
   west = - (los - w*cos(ze))/sin(ze)
   westTemp  = fpiT(l)
   westTflag = fpiTflag(l)
   westVflag = fpiVflag(l)
endif else begin
   west = [0.0]
   westTime = [min(fpitime(0))]
   westTemp = [0.0]
   westTflag = [0.0]
   westVflag = [0.0]
endelse

fpi = {NorthTime:NorthTime, $
       SouthTime:SouthTime, $
       EastTime:EastTime, $
       WestTime:WestTime, $
       ZenithTime:ZenithTime, $
       north: north, $
       south:south, $
       east:east, $
       west:west, $
       northTemp: northTemp, $
       southTemp:southTemp, $
       eastTemp:eastTemp, $
       westTemp:westTemp, $
       northTflag: northTflag, $
       southTflag:southTflag, $
       eastTflag:eastTflag, $
       westTflag:westTflag, $
       northVflag: northVflag, $
       southVflag:southVflag, $
       eastVflag:eastVflag, $
       westVflag:westVflag, $
       fpitime:fpitime}

end


