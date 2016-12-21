
pro read_cdf_imf, magfile, swfile, magfield, sw, magtime, swtime, isgse, $
                  delay, xyz

  ;---------------------------------------------------------------------------
  ; mag file first.
  ;---------------------------------------------------------------------------

  id = CDF_OPEN(magfile)
  
  info = cdf_inquire(id)
  
  natt = info.natts
  nvar = info.nvars
  if nvar eq 0 then begin
    nvar = info.nzvars
    zvar = 1
  endif else zvar = 0
  
  vars = strarr(nvar)
  
  nrows = info.maxrec
  
  if (nrows lt 1) then begin
    cdf_control, id, variable=0, /zvariable,get_var_info=v
    nrows = v.maxrec
  endif
  
  xyz      = fltarr(3, nrows) - 1.0e32
  magfield = fltarr(3, nrows) - 1.0e32
  magtime  = dblarr(nrows)

  epoch =  -1
  imf   =  -1
  pos   =  -1
  isgse =   0
  delay = 0.0

  for i=0,nvar-1 do begin
  
    if (zvar eq 0) then r = cdf_varinq(id,i)		$
    else r = cdf_varinq(id,i,/zvariable)
    vars(i) = r.name
    tmp = mklower(vars(i))
    if (strpos(tmp,'epoch') gt -1 and epoch lt 0) then epoch = i
    if (strpos(tmp,'p') gt -1 and $
        strpos(tmp,'gsm') gt -1 and $
        pos eq -1) then pos = i
    if (strpos(tmp,'x-gsm') gt -1 and $
        pos eq -1) then pos = i
    if (strpos(tmp,'b') gt -1 and $
        strpos(tmp,'gsm') gt -1  and $
        imf eq -1) then imf = i

    if (strpos(tmp,'b') gt -1 and $
        strpos(tmp,'gsm') gt -1 and $
        strpos(tmp,'label') lt 0 and imf eq -1) then imf = i

    ;print, tmp, epoch, pos, imf

  endfor

  swap = 0
  if (epoch gt -1) then begin
    cdf_varget, id, vars(epoch), times, rec_count=nrows
    if (max(times) lt 1.0e10) then swap = 1
    print, mm(times), swap
    magtime = reform(times(0,*))
  endif

  if (pos gt -1) then begin
    cdf_varget, id, vars(pos), xyz, rec_count=nrows
    PosInMagFile = 1
    ; pos = -1
  endif else begin
    print, 'Could not find position array within the magfile.'
    PosInMagFile = 0
  endelse

  if (imf gt -1) then begin
      tmp = mklower(vars(imf))
      if (strpos(tmp,'gse') gt -1) then isgse = 1
      cdf_varget, id, vars(imf), magfield, rec_count=nrows
  endif else begin
      print, 'Could not find position array. Please check CDF file!'
      print, 'Here are the variables within this file:'
      for i=0,nvar-1 do print, vars(i)
  endelse

  if (swap) then begin
      magtime  = swap_endian(magtime)
      xyz      = swap_endian(xyz)
      magfield = swap_endian(magfield)
  endif

  cdf_close, id 

  ;---------------------------------------------------------------------------
  ; solar wind file next.
  ;---------------------------------------------------------------------------

  if (strpos(swfile,"cdf") lt 0) then begin

      openr,1,swfile

      line = ""

      readf,1,line
      readf,1,line

      if (strpos(line, 'GLOBAL') gt 0) then begin

          while (strpos(line,'dd-mm-yyyy') lt 0) do readf,1,line
          i = 0
          while not eof(1) do begin
              readf,1,line
              i = i + 1
          endwhile
          nLines = i - 2

          close,1
          openr,1,swfile
          while (strpos(line,'dd-mm-yyyy') lt 0) do readf,1,line

          dd = 0
          mm = 0
          yy = 0
          hh = 0
          mi = 0
          ss = 0
          np = 0.0
          pr = 0.0
          vx = 0.0
          vy = 0.0
          vz = 0.0

          swtime = dblarr(nLines)
          sw     = fltarr(6,nLines)

          for iLine = 0, nLines-1 do begin

              readf,1,line
              yy = strmid(line,6,4)
              mm = strmid(line,3,2)
              dd = strmid(line,0,2)
              hh = strmid(line,11,2)
              mi = strmid(line,14,2)
              ss = strmid(line,17,2)
              np = strmid(line,29,8)
              pr = strmid(line,37,14)
              vx = strmid(line,51,14)
              vy = strmid(line,65,14)
              vz = strmid(line,79,14)
              
              cdf_epoch, rtime, yy,mm,dd,hh,mi,ss,/compute

              swtime(iLine) = rtime
              sw(0:2,iLine) = [vx, vy, vz]
              sw(3,iLine) = np
              sw(4,iLine) = 100000.0
              sw(5,iLine) = 0.1
              print, yy,mm,dd,hh,mi,ss, iLine, nLines, vx, vy, vz, sw(4,iLine)

          endfor

          close,1

      endif else begin

          for i=0,41 do readf,1,line
          iyear = 0
          day = 0d
          np = 0.0
          tp = 0.0
          speed = 0.0
          alpharatio = 0.0
          vx = 0.0
          vy = 0.0
          vz = 0.0
          itime = intarr(6)
          swtime = dblarr(648)
          sw     = fltarr(5,648)
          for i=0,647 do begin
              readf,1,iyear, day, np, tp, speed, alpharatio, $
                vx, vy, vz
              
              itime(0) = iyear
              itime(1) = 1
              itime(2) = fix(day)
              nh = (day - itime(2)) * 24.0d
              itime(3) = fix(nh)
              nm = (nh - itime(3))*60.0d
              itime(4) = fix(nm)
              ns = (nm - itime(4))*60.0d
              itime(5) = fix(ns)

              day2 = (day-fix(day))*86400.0
              
              c_a_to_r, itime, rtime
              c_r_to_a, itime, rtime

              cdf_epoch, rtime, $
                itime(0), itime(1), itime(2), itime(3), itime(4), itime(5),/compute

              swtime(i) = rtime

              if (vx lt -5000.0) then begin
                  if (speed lt 0.0) then begin
                      vx = sw(0,i-1)
                      vy = sw(1,i-1)
                      vz = sw(2,i-1)
                  endif else begin
                      vxold = sw(0,i-1)
                      vyold = sw(1,i-1)
                      vzold = sw(2,i-1)
                      speedold = sqrt(vxold^2+vyold^2+vzold^2)
                      vx = speed * vxold/speedold
                      vy = speed * vyold/speedold
                      vz = speed * vzold/speedold
                  endelse
              endif
              if (np lt 0.0) then np = sw(3,i-1)
              if (tp lt 0.0) then np = sw(4,i-1)
              sw(0:2,i) = [vx, vy, vz]
              sw(3,i) = np
              sw(4,i) = tp
          endfor
          close,1

      endelse

  endif else begin

      print, "solar wind file : ", swfile

      id = CDF_OPEN(swfile)
  
      info = cdf_inquire(id)
  
      natt = info.natts
      nvar = info.nvars
      if nvar eq 0 then begin
          nvar = info.nzvars
          zvar = 1
      endif else zvar = 0
  
      vars = strarr(nvar)
  
      nrows = info.maxrec
  
      if (nrows lt 1) then begin
          cdf_control, id, variable=0, /zvariable,get_var_info=v
          nrows = v.maxrec
      endif
  
      sw     = fltarr(6, nrows) - 1.0e32
      swtime = dblarr(nrows)
      epoch  =  -1
      Np     =  -1
      Vp     =  -1
      Tp     =  -1
      iAlpha =  -1
      delay  = 0.0
      swmag  =   0

      for i=0,nvar-1 do begin
  
          if (zvar eq 0) then r = cdf_varinq(id,i)		$
          else r = cdf_varinq(id,i,/zvariable)
          vars(i) = r.name
          tmp = mklower(vars(i))
          if (strpos(tmp,'epoch') gt -1) then epoch = i
          if (strpos(tmp,'p') gt -1 and $
              strpos(tmp,'gs') gt -1 and $
              pos eq -1) then pos = i

;          if (((strpos(tmp,'v') gt -1 and $
;                strpos(tmp,'gs') gt -1) or $ ; for WIND and IMP-8
;               (tmp eq 'vp')) and $ ; for ACE
;              Vp eq -1) then begin
;              Vp = i
;              if (tmp eq 'vp') then swmag = 1
;          endif

          if ((strpos(tmp,'v') gt -1 and $
               strpos(tmp,'gsm') gt -1) and $ 
              Vp eq -1) then begin
              Vp = i
          endif

          if (strpos(tmp,'el_bulk_vel_magn') gt -1 and $ ; for wind 
              Vp eq -1) then begin
              Vp = i
          endif

          if (((strpos(tmp,'thermal') gt -1) or $ ; for WIND and IMP-8
               (tmp eq 'tpr') or (tmp eq 'te')) and $ ; for ACE
              Tp eq -1) then Tp = i

          if (tmp eq 'np' and $ ; for ACE and WIND and IMP-8
              Np eq -1) then Np = i

          if (tmp eq 'el_density' and $ ; for WIND
              Np eq -1) then Np = i

          if (strpos(tmp,'alpha') gt -1 and $ ; for ACE and WIND and IMP-8
              iAlpha eq -1) then iAlpha = i

;          print, vars(i), epoch, pos, vp, tp, np, iAlpha

      endfor

      if (epoch gt -1) then begin
          cdf_varget, id, vars(epoch), times, rec_count=nrows
          swtime = reform(times(0,*))
          if (swap) then swtime = swap_endian(swtime)
      endif

      if (pos gt -1 and not PosInMagFile) then begin
          cdf_varget, id, vars(pos), xyz, rec_count=nrows
          if (swap) then xyz = swap_endian(xyz)
      endif else begin
          if (pos eq -1) then begin
              print, 'Could not find position array within the swfile.'
              print, 'Can not calculate the time delay.'
          endif
      endelse

      if (Vp gt -1) then begin
          tmp = mklower(vars(Vp))
          cdf_varget, id, vars(Vp), velocity, rec_count=nrows
          if (swap) then velocity = swap_endian(velocity)
          ;print, "v:",mm(velocity)

          ;help, velocity
          nX = n_elements(velocity(*,0))

          if (nX eq 3) then $  ; This means that the satellite is probably ACE
             sw(0:2,*) = velocity $
          else begin
             sw(0,*) = -abs(velocity)
             sw(1,*) = 0.0
             sw(2,*) = 0.0
             endelse                    
      endif else begin
          print, 'Could not find velocity array. Please check CDF file!'
          print, 'Here are the variables within this file:'
          for i=0,nvar-1 do print, vars(i)
      endelse

      if (Np gt -1) then begin
          tmp = mklower(vars(Np))
          cdf_varget, id, vars(Np), numden, rec_count=nrows
          if (swap) then numden = swap_endian(numden)
          sw(3,*) = numden(0,*)
      endif else begin
          print, 'Could not find number density array. Please check CDF file!'
          print, 'Here are the variables within this file:'
          for i=0,nvar-1 do print, vars(i)
      endelse

      if (iAlpha gt -1) then begin
          tmp = mklower(vars(iAlpha))
          cdf_varget, id, vars(iAlpha), alpha, rec_count=nrows
          if (swap) then alpha = swap_endian(alpha)
          sw(5,*) = alpha(0,*)
      endif else begin
          print, 'Could not find alpha ratio array. Please check CDF file!'
          print, 'Here are the variables within this file:'
          for i=0,nvar-1 do print, vars(i)
      endelse

      if (Tp gt -1) then begin
          tmp = mklower(vars(Tp))
          cdf_varget, id, vars(Tp), temperature, rec_count=nrows

          if (swap) then temperature = swap_endian(temperature)

          for i=1,nrows-1 do if temperature(0,i) lt 0.0 then $
            temperature(0,i) = temperature(0,i-1)

          sw(4,*) = temperature(0,*)

          goodsw = where(abs(sw(4,*)) lt 2000.0, count)
          if (count gt 0) then begin
              if (mean(sw(4,goodsw)) lt 100.0) then begin
                  mik = 1.6726e-27/1.3807e-23
                  sw(4,goodsw) = ((sw(4,goodsw)*1000.0)^2.0)*mik
              endif
          endif

      endif else begin
          print, 'Could not find temperature array. Please check CDF file!'
          print, 'Here are the variables within this file:'
          for i=0,nvar-1 do print, vars(i)
      endelse

      cdf_close, id 

      l = where(reform(sw(4,*)) gt 0.0 and reform(sw(0,*)) gt -1500, c)
      if (c gt 0) then begin
          nrows = c
          swnew = fltarr(6,c)
          for i=0,5 do swnew(i,*) = sw(i,l)
          sw = swnew
          swtime = swtime(l)
      endif

  endelse

  if (abs(xyz(0,0)) lt 1.0e31 and min(abs(sw(0,*))) lt 1.0e31) then begin
    good = where(abs(xyz(0,*)) lt 300.0*6372.0)
    r = sqrt(xyz(0,good)^2 + xyz(1,good)^2 + xyz(2,good)^2)
    loc = where(r lt 6372.0, count)
    if (count gt 0) then xyz(*,good) = xyz(*,good)*6372.0
    goodsw = where(abs(sw(0,*)) lt 2000.0, count)
    if (count gt 0) then begin
      delay = (mean(xyz(0,good))-32.0*6372.0) / mean(abs(sw(0,goodsw)))
    endif else delay = (mean(xyz(0,good))-32.0*6372.0) / 400.0
    if (delay gt 7200.0 or delay lt -1800.0) then begin
      print, 'Time delay looks wrong : ',delay
      print, 'Setting to 0'
      delay = 0.0
    endif
  endif

end


