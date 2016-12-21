
pro get_imf, stime, etime, data, dt = dt

  if (n_elements(dt) eq 0) then dt = -1

  if (n_elements(stime) eq 6) then begin
    itime = stime
    c_a_to_r, itime, stime
  endif

  if (n_elements(etime) eq 6) then begin
    itime = etime
    c_a_to_r, itime, etime
  endif

  cTime = sTime - 7200.0

  IsFirstTime = 1

  while (cTime lt eTime+24.0*3600.0-1.0) do begin

    c_r_to_a, iTime, cTime
    syear  = tostr(iTime(0))
    smonth = chopr('0'+tostr(iTime(1)),2)
    sday   = chopr('0'+tostr(iTime(2)),2)

    dir    = getenv('ACE_DATA')+'/'+syear+'/'
    file   = dir+'ace'+syear+smonth+sday+'.save'

    print, 'restoring :',file
    restore, file
    delay = 3600.0

    if (IsFirstTime) then begin
      bSave = b
      vSave = v
      bmagSave = reform(bmag)
      eklSave  = reform(ekl)
      maSave   = reform(ma)
      vaSave   = reform(va)
      denSave  = n
      pbSave   = reform(pb)
      pdynSave = pdyn
      tSave    = t
      timeSave = time + delay
      IsFirstTime = 0
    endif else begin

      nsave = n_elements(bSave(0,*))
      nold  = n_elements(b(0,*))

      bSaveNew = fltarr(3,nsave+nold)
      bSaveNew(*,0:nsave-1) = bSave
      bSaveNew(*,nsave:nsave+nold-1) = b
      bSave = bSaveNew

      nsave = n_elements(vSave(0,*))
      nold  = n_elements(v(0,*))

      vSaveNew = fltarr(3,nsave+nold)
      vSaveNew(*,0:nsave-1) = vSave
      vSaveNew(*,nsave:nsave+nold-1) = v
      vSave = vSaveNew

      bmagSave = [bmagSave,reform(bmag)]
      eklSave  = [eklSave,reform(ekl)]
      maSave   = [maSave,reform(ma)]
      vaSave   = [vaSave,reform(va)]
      denSave  = [denSave, n]
      pbSave   = [pbSave,reform(pb)]
      pdynSave = [pdynSave,pdyn]
      tSave    = [tSave, t]
      timeSave = [timeSave,time]

    endelse

    cTime = cTime + 24.0*3600.0

  endwhile

  loc = where(timesave ge sTime and timesave lt eTime, count)

  if (count gt 0) then begin
    data = { bx   : reform(bSave(0,loc)), $
             by   : reform(bSave(1,loc)), $
             bz   : reform(bSave(2,loc)), $
             bmag : bmagSave(loc),        $
             vx   : reform(vSave(0,loc)), $
             vy   : reform(vSave(1,loc)), $
             vz   : reform(vSave(2,loc)), $
             Ekl  : eklSave(loc),         $
             mach : maSave(loc),          $
             va   : vaSave(loc),          $
             den  : denSave(loc),         $
             DynP : pdynSave(loc),        $
             MagP : pbSave(loc),          $
             Temp : tSave(loc),           $
             Time : timeSave(loc)}

    if (dt gt 0.0) then begin

      nTimes = (etime-stime)/dt + 1
      datanew = { bx   :  fltarr(nTimes), $
                  by   :  fltarr(nTimes), $
             	  bz   :  fltarr(nTimes), $
             	  bmag :  fltarr(nTimes), $
             	  vx   :  fltarr(nTimes), $
             	  vy   :  fltarr(nTimes), $
             	  vz   :  fltarr(nTimes), $
             	  Ekl  :  fltarr(nTimes), $
             	  mach :  fltarr(nTimes), $
             	  va   :  fltarr(nTimes), $
             	  den  :  fltarr(nTimes), $
             	  DynP :  fltarr(nTimes), $
             	  MagP :  fltarr(nTimes), $
             	  Temp :  fltarr(nTimes), $
             	  Time :  dblarr(nTimes)}

      print, "Interpolating IMF onto time grid..."

      Nearest = 1
      if (Nearest) then begin

        dn = nTimes/10
        iStart = 0
        iEnd   = n_elements(data.time)-1

        l = 1L

        for i=0L,nTimes-1 do begin

          if (fix(i/dn) ne fix((i-1)/dn)) then print, float(i*10)/dn,"% done"
          t = sTime + i*dt

          while (data.time(l) lt t and l lt iEnd) do l = l + 1L

          x = 1.0 - (data.time(l) - t)/ (data.time(l) - data.time(l-1))

  	  datanew.bx(i)   = x*data.bx(l) + (1-x)*data.bx(l-1)
  	  datanew.by(i)   = x*data.by(l) + (1-x)*data.by(l-1)
  	  datanew.bz(i)   = x*data.bz(l) + (1-x)*data.bz(l-1)
  	  datanew.bmag(i) = x*data.bmag(l) + (1-x)*data.bmag(l-1)
  	  datanew.vz(i)   = x*data.vx(l) + (1-x)*data.vx(l-1)
  	  datanew.vy(i)   = x*data.vy(l) + (1-x)*data.vy(l-1)
  	  datanew.vz(i)   = x*data.vz(l) + (1-x)*data.vz(l-1)
  	  datanew.Ekl(i)  = x*data.Ekl(l) + (1-x)*data.ekl(l-1)
  	  datanew.mach(i) = x*data.mach(l) + (1-x)*data.mach(l-1)
  	  datanew.va(i)   = x*data.va(l) + (1-x)*data.va(l-1)
  	  datanew.den(i)  = x*data.den(l) + (1-x)*data.den(l-1)
  	  datanew.DynP(i) = x*data.DynP(l) + (1-x)*data.DynP(l-1)
  	  datanew.MagP(i) = x*data.MagP(l) + (1-x)*data.MagP(l-1)
  	  datanew.Temp(i) = x*data.Temp(l) + (1-x)*data.Temp(l-1)
          datanew.time(i) = t
        endfor
      endif else begin

        for i=0L,nTimes-1 do begin
          t = sTime + i*dt
          loc = where(data.Time ge t-dt/2.0 and data.Time lt t+dt/2.0,count)
          if (count eq 0) then begin
            d = abs(data.Time - t)
            loc = where(d eq min(d))
          endif
          datanew.bx(i)   = mean(data.bx(loc))
          datanew.by(i)   = mean(data.by(loc))
          datanew.bz(i)   = mean(data.bz(loc))
          datanew.bmag(i) = mean(data.bmag(loc))
          datanew.vz(i)   = mean(data.vx(loc))
          datanew.vy(i)   = mean(data.vy(loc))
          datanew.vz(i)   = mean(data.vz(loc))
          datanew.Ekl(i)  = mean(data.Ekl(loc))
          datanew.mach(i) = mean(data.mach(loc))
          datanew.va(i)   = mean(data.va(loc))
          datanew.den(i)  = mean(data.den(loc))
          datanew.DynP(i) = mean(data.DynP(loc))
          datanew.MagP(i) = mean(data.MagP(loc))
          datanew.Temp(i) = mean(data.Temp(loc))
        endfor
      endelse

      data = datanew

    endif

  endif else begin
    data = { bx   : 0.0, $
             by   : 0.0, $
             bz   : 0.0, $
             bmag : 0.0, $
             vx   : 0.0, $
             vy   : 0.0, $
             vz   : 0.0, $
             Ekl  : 0.0, $
             mach : 0.0, $
             va   : 0.0, $
             den  : 0.0, $
             DynP : 0.0, $
             MagP : 0.0, $
             Temp : 0.0, $
             Time : 0.0}
  endelse


end
